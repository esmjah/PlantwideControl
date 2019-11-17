# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 17:00:14 2014

@author: codas
"""

import casadi as C
import numpy as NP
import scipy
import scipy.io as sio
from scipy.interpolate import interp1d




class DaeModel:
    
    
    ocp = None
    DT = None
    measurmentList = None
    mSXF = None
    G = None 
    z0 = None
    
    measurementScaling = None
    stateScaling = None
    algStateScaling = None
    

    def init(self,ocp,DT,measuremntsList):
        
        measurementScaling = C.vertcat([ocp.variable(k).nominal for k in measuremntsList])
        stateScaling = C.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])   
        algStateScaling = C.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())]) 
        controlScaling = C.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())]) 
        
        odeS = C.substitute(ocp.ode(ocp.x),C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/stateScaling
        algS = C.substitute(ocp.alg,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))
        
        sysIn = C.daeIn(x=ocp.x,z=ocp.z,p=ocp.u,t=ocp.t)
        sysOut = C.daeOut(ode=odeS,alg=algS)
        odeF=C.SXFunction(sysIn,sysOut)    
        odeF.init()
        
        
        C.Integrator.loadPlugin("idas")
        G = C.Integrator("idas",odeF)
        G.setOption("reltol",1e-8)  #for IDAS
        G.setOption("abstol",1e-8)  #for IDAS 
        G.setOption("max_multistep_order",5)  #for IDAS
        G.setOption("tf",DT)		        			
        G.init()

        mSX = C.vertcat([ocp.variable(measuremntsList[k]).beq  for k in range(len(measuremntsList))])
        mSX = C.substitute(mSX,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/measurementScaling
        mSXF = C.SXFunction(sysIn,[mSX])
        mSXF.init()
        
        self.measurementScaling = measurementScaling
        self.stateScaling = stateScaling
        self.algStateScaling = algStateScaling
        self.controlScaling = controlScaling


        self.ocp = ocp
        self.DT = DT
        self.measurmentList = measuremntsList
        self.mSXF = mSXF
        self.G = G
        
        self.z0 = C.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])/algStateScaling
        self.x0 = C.vertcat([ocp.variable(ocp.x[k].getName()).initialGuess.getValue() for k in range(ocp.x.size())])/stateScaling
        
    def computeJacobians(self,x,u,measurmentTags): 

  
        ocp = self.ocp        
        
        sysIn = C.daeIn(x=ocp.x,z=ocp.z,p=ocp.u,t=ocp.t)
        sysOut = C.daeOut(ode=ocp.ode(ocp.x),alg=ocp.alg)
        odeF=C.SXFunction(sysIn,sysOut)    
        odeF.init()

        mSX = C.vertcat([ocp.variable(measurmentTags[k]).beq  for k in range(len(measurmentTags))])
        mSXF = C.SXFunction(sysIn,[mSX])
        mSXF.init()
     
        AS = odeF.jac('x','ode')
        BS = odeF.jac('p','ode')
        CS = mSXF.jac('x',0)
        DS = mSXF.jac('p',0)
        
        
        funcJacs = C.SXFunction(sysIn,[AS,BS,CS,DS])
        funcJacs.init()
        
              
        funcJacs.setInput(x,'x')
        funcJacs.setInput(u,'p')
        
        funcJacs.evaluate()
        
        
        a = funcJacs.getOutput(0)
        b = funcJacs.getOutput(1)
        c = funcJacs.getOutput(2)
        d = funcJacs.getOutput(3)        
        
        return a,b,c,d
        
    def checkStability(self,x,u):
        
        stateScaling =  self.stateScaling
        algStateScaling = self.algStateScaling
        controlScaling = self.controlScaling
        ocp = self.ocp        
        
        odeS = C.substitute(ocp.ode(ocp.x),C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/stateScaling
        algS = C.substitute(ocp.alg,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))
        
        sysIn = C.daeIn(x=ocp.x,z=ocp.z,p=ocp.u,t=ocp.t)
        sysOut = C.daeOut(ode=odeS,alg=algS)
        odeF=C.SXFunction(sysIn,sysOut)    
        odeF.init()
        
        dG = odeF.jacobian('x','ode')
        dG.init()
        dG.setInput(x/stateScaling,'x')
        dG.setInput(u/controlScaling,'p')
        
        dG.evaluate()
        
        dGdx = dG.getOutput()
        
        eigVals, eigVecs = scipy.linalg.eig(dGdx)
        
        return eigVals,dGdx,odeS
        
        
    def changeMeasurement(self,measuremntsList):
        measurementScaling = C.vertcat([self.ocp.variable(k).nominal for k in measuremntsList])
        ocp = self.ocp
        mSX = C.vertcat([ocp.variable(measuremntsList[k]).beq  for k in range(len(measuremntsList))])
        mSX = C.substitute(mSX,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([self.stateScaling*ocp.x,self.algStateScaling*ocp.z,self.controlScaling*ocp.u]))/measurementScaling
        sysIn = C.daeIn(x=ocp.x,z=ocp.z,p=ocp.u,t=ocp.t)
        
        mSXF = C.SXFunction(sysIn,[mSX])
        mSXF.init()

        self.measurementScaling = measurementScaling
        self.mSXF = mSXF


    def oneStep(self,x0,u,z0=None):
        
        if self.ocp.z.size() > 0:
            if z0 is None:        
                z0 = self.z0
            else:
                z0 = z0/self.algStateScaling
        else:
            z0 = self.z0
            
        x0 = x0/self.stateScaling
        u = u/self.controlScaling

                     
        
        self.G.setInput(x0,'x0')
        self.G.setInput(z0,'z0')
        self.G.setInput(u,'p')
            
        self.G.evaluate()
            
        xk = self.G.getOutput('xf')
        zf = self.G.getOutput('zf')
        
        
        self.mSXF.setInput(xk,'x')
        self.mSXF.setInput(zf,'z')
        self.mSXF.setInput(u,'p')
        
        self.mSXF.evaluate()
        
        y = self.mSXF.getOutput()*self.measurementScaling
        xk = xk*self.stateScaling
        zf = zf*self.algStateScaling
                        
        return xk,zf,y  
        
    def findZ(self,u,x0=None,z0=None):
        
        if z0 is None:        
            z0 = self.z0
        else:
            z0 = z0/self.algStateScaling
        if x0 is None:    
            x0 = self.x0
        else:
            x0 = x0/self.stateScaling
            
        u = u/self.controlScaling
              
        ocp = self.ocp
        
        stateScaling = self.stateScaling   
        algStateScaling = self.algStateScaling
        controlScaling = self.controlScaling 
         
        
        algS = C.substitute(ocp.alg,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))
                
        nlp=C.SXFunction(C.nlpIn(x=ocp.z,p=C.vertcat([ocp.x,ocp.u])),C.nlpOut(f=C.SX(0),g=algS))
        C.NlpSolver.loadPlugin('ipopt')
        solver = C.NlpSolver('ipopt',nlp)
        solver.setOption('print_user_options','no')
        solver.setOption('tol',1e-10)
        solver.setOption('print_level',0)
        solver.setOption('file_print_level',0)
        solver.setOption("max_iter",200)	# IPOPT maximum iterations

      
        solver.init()
        
        solver.setInput(NP.zeros(ocp.z.size()),'lbg')	# g_L
        solver.setInput(NP.zeros(ocp.z.size()),'ubg')	# g_U        
  

        solver.setInput(z0,'x0')   
        solver.setInput(C.vertcat([x0,u]), 'p')       
  
        solver.evaluate()

        self.z0 = solver.output('x')
        
        return solver.output('x')*self.algStateScaling
        


    def findSteadyState(self,u,x0,z0=None,simCount=0,consList=[]):
        
        if z0 is not None:        
            z0 = z0/self.algStateScaling
        else:
            z0 = self.z0
        if x0 is not None:    
            x0 = x0/self.stateScaling
        else:
            x0 = self.x0
                
        ocp = self.ocp
        u = u/self.controlScaling

        measurementScaling = self.measurementScaling
        stateScaling = self.stateScaling   
        algStateScaling = self.algStateScaling 
        controlScaling = self.controlScaling 
        
        consScaling = C.vertcat([ocp.variable(k).nominal for k in consList])
        
  
        for k in range(simCount):
            x0,z0,y = self.oneStep(x0*stateScaling,u*controlScaling,z0*algStateScaling)
            
            x0 = x0/stateScaling
            z0 = z0/algStateScaling
            y = y/measurementScaling
        
        
        odeS = C.substitute(ocp.ode(ocp.x),C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/stateScaling
        algS = C.substitute(ocp.alg,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))
                                      
  
        mSX = C.vertcat([ocp.variable(consList[k]).beq  for k in range(len(consList))])
        mSX = C.substitute(mSX,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/consScaling
  
  
        split = C.SXFunction([C.vertcat([ocp.x,ocp.z])],[ocp.x,ocp.z])
        split.init()
      
        nlp=C.SXFunction(C.nlpIn(x=C.vertcat([ocp.x,ocp.z]),p=ocp.u),C.nlpOut(f=C.SX(0),g=C.vertcat([odeS,algS,mSX])))
        C.NlpSolver.loadPlugin('ipopt')
        solver = C.NlpSolver('ipopt',nlp)
        solver.setOption('print_user_options','no')
        solver.setOption('tol',1e-10)
        solver.setOption('print_level',5)
        solver.setOption("max_iter",2000)	# IPOPT maximum iterations
              
        solver.init()

        xMin = C.vertcat([ocp.variable(ocp.x[i].getName()).min.getValue() for i in range(ocp.x.size())])/stateScaling
        xMax = C.vertcat([ocp.variable(ocp.x[i].getName()).max.getValue() for i in range(ocp.x.size())])/stateScaling
     
        zMin = C.vertcat([ocp.variable(ocp.z[i].getName()).min.getValue() for i in range(ocp.z.size())])/algStateScaling
        zMax = C.vertcat([ocp.variable(ocp.z[i].getName()).max.getValue() for i in range(ocp.z.size())])/algStateScaling
        
        cMin = C.vertcat([ocp.variable(consList[i]).min.getValue() for i in range(len(consList))])/consScaling
        cMax = C.vertcat([ocp.variable(consList[i]).max.getValue() for i in range(len(consList))])/consScaling 
        
        solver.setInput(C.vertcat([NP.zeros(ocp.z.size()+ocp.x.size()),cMin]),'lbg')	# g_L
        solver.setInput(C.vertcat([NP.zeros(ocp.z.size()+ocp.x.size()),cMax]),'ubg')	# g_U        
  
  
        solver.setInput(C.vertcat([xMin,zMin]),  'lbx')					# u_L
        solver.setInput(C.vertcat([xMax,zMax]),  'ubx')					# u_U 
  

        solver.setInput(C.vertcat([x0,z0]),'x0')   
        solver.setInput(u, 'p')       
  
        solver.evaluate()


        xz = solver.output('x')
        split.setInput(xz)
        split.evaluate()
        x0 = split.getOutput(0)
        z0 = split.getOutput(1)
        
        self.mSXF.setInput(x0,'x')
        self.mSXF.setInput(z0,'z')
        self.mSXF.setInput(u,'p')
        
        self.mSXF.evaluate()
        
        y0 = self.mSXF.getOutput()        
        
        
        return x0*stateScaling,z0*algStateScaling,y0*measurementScaling




    def findOptimalPoint(self,u0,x0=None,z0=None,simCount=0,consList=[]):
        
        if z0 is not None:        
            z0 = z0/self.algStateScaling
        else:
            z0 = self.z0
        if x0 is not None:    
            x0 = x0/self.stateScaling
        else:
            x0 = self.x0
                
        ocp = self.ocp
        u0 = u0/self.controlScaling

        measurementScaling = self.measurementScaling
        stateScaling = self.stateScaling   
        algStateScaling = self.algStateScaling 
        controlScaling = self.controlScaling 
        
        consScaling = C.vertcat([ocp.variable(k).nominal for k in consList])
        
  
        for k in range(simCount):
            x0,z0,y = self.oneStep(x0*stateScaling,u*controlScaling,z0*algStateScaling)
            
            x0 = x0/stateScaling
            z0 = z0/algStateScaling
            y = y/measurementScaling
        
        
        odeS = C.substitute(ocp.ode(ocp.x),C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/stateScaling
        algS = C.substitute(ocp.alg,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))
        ltermS = C.substitute(ocp.lterm,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))
                                     
  
        mSX = C.vertcat([ocp.variable(consList[k]).beq  for k in range(len(consList))])
        mSX = C.substitute(mSX,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/consScaling
  
  
        pathVarNames = [sv.getName()  for sv in ocp.beq(ocp.path)]
        pathScaling = C.vertcat([ocp.nominal(pv)  for pv in pathVarNames])        
        pathS = C.substitute(ocp.beq(ocp.beq(ocp.path)),C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/pathScaling
  
  
        split = C.SXFunction([C.vertcat([ocp.u,ocp.x,ocp.z])],[ocp.u,ocp.x,ocp.z])
        split.init()
      
        nlp=C.SXFunction(C.nlpIn(x=C.vertcat([ocp.u,ocp.x,ocp.z])),C.nlpOut(f=ltermS,g=C.vertcat([odeS,algS,mSX,pathS])))
        C.NlpSolver.loadPlugin('ipopt')
        solver = C.NlpSolver('ipopt',nlp)
        solver.setOption('print_user_options','no')
        solver.setOption('print_level',5)
        solver.setOption("tol",1e-14) 
        solver.setOption("max_iter",300)	# IPOPT maximum iterations
              
        solver.init()
        
        uMin = C.vertcat([ocp.variable(ocp.u[i].getName()).min.getValue() for i in range(ocp.u.size())])/controlScaling
        uMax = C.vertcat([ocp.variable(ocp.u[i].getName()).max.getValue() for i in range(ocp.u.size())])/controlScaling        

        xMin = C.vertcat([ocp.variable(ocp.x[i].getName()).min.getValue() for i in range(ocp.x.size())])/stateScaling
        xMax = C.vertcat([ocp.variable(ocp.x[i].getName()).max.getValue() for i in range(ocp.x.size())])/stateScaling
     
        zMin = C.vertcat([ocp.variable(ocp.z[i].getName()).min.getValue() for i in range(ocp.z.size())])/algStateScaling
        zMax = C.vertcat([ocp.variable(ocp.z[i].getName()).max.getValue() for i in range(ocp.z.size())])/algStateScaling
        
        cMin = C.vertcat([ocp.variable(consList[i]).min.getValue() for i in range(len(consList))])/consScaling
        cMax = C.vertcat([ocp.variable(consList[i]).max.getValue() for i in range(len(consList))])/consScaling 
        
           
        
        pathMax = C.vertcat([ocp.variable(pathVarNames[i]).max.getValue()  for i in range(ocp.path.size())])/pathScaling;
        pathMin = C.vertcat([ocp.variable(pathVarNames[i]).min.getValue()  for i in range(ocp.path.size())])/pathScaling;
    
        
        
        
        solver.setInput(C.vertcat([NP.zeros(ocp.z.size()+ocp.x.size()),cMin,pathMin]),'lbg')	# g_L
        solver.setInput(C.vertcat([NP.zeros(ocp.z.size()+ocp.x.size()),cMax,pathMax]),'ubg')	# g_U        
  
  
        solver.setInput(C.vertcat([uMin,xMin,zMin]),  'lbx')					# u_L
        solver.setInput(C.vertcat([uMax,xMax,zMax]),  'ubx')					# u_U 
  

        solver.setInput(C.vertcat([u0,x0,z0]),'x0')   
  
        solver.evaluate()


        xz = solver.output('x')
        split.setInput(xz)
        split.evaluate()
        u0 = split.getOutput(0)
        x0 = split.getOutput(1)
        z0 = split.getOutput(2)
        
        self.mSXF.setInput(x0,'x')
        self.mSXF.setInput(z0,'z')
        self.mSXF.setInput(u0,'p')
        
        self.mSXF.evaluate()
        
        y0 = self.mSXF.getOutput()        
        
        
        return u0*controlScaling,x0*stateScaling,z0*algStateScaling,y0*measurementScaling

















        
    def findConstrainedSteadyState(self,altUTags,altUVals,u0,x0=None,z0=None,consList=[]):
        
        if z0 is not None:        
            z0 = z0/self.algStateScaling
        else:
            z0 = self.z0
        if x0 is not None:    
            x0 = x0/self.stateScaling
        else:
            x0 = self.x0
    
        u0 = u0/self.controlScaling

                
        ocp = self.ocp

        measurementScaling = self.measurementScaling
        stateScaling = self.stateScaling   
        algStateScaling = self.algStateScaling 
        controlScaling = self.controlScaling 
        
        consScaling = C.vertcat([ocp.variable(k).nominal for k in consList])
        altUScaling = C.vertcat([ocp.variable(altUTags[k]).nominal for k in range(len(altUTags))]) 
        
        odeS = C.substitute(ocp.ode(ocp.x),C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/stateScaling
        algS = C.substitute(ocp.alg,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))
                                      
        altU = C.vertcat([ocp.variable(altUTags[k]).beq  for k in range(len(altUTags))])
        altU = C.substitute(altU,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/altUScaling
    
        mSX = C.vertcat([ocp.variable(consList[k]).beq  for k in range(len(consList))])
        mSX = C.substitute(mSX,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/consScaling
  
  
        split = C.SXFunction([C.vertcat([ocp.x,ocp.z,ocp.u])],[ocp.x,ocp.z,ocp.u])
        split.init()
      
        nlp=C.SXFunction(C.nlpIn(x=C.vertcat([ocp.x,ocp.z,ocp.u])),C.nlpOut(f=C.SX(0),g=C.vertcat([odeS,algS,altU,mSX])))
        C.NlpSolver.loadPlugin('ipopt')
        solver = C.NlpSolver('ipopt',nlp)
        solver.setOption('tol',1e-10)
        solver.setOption('print_user_options','yes')
        solver.setOption("max_iter",200)	# IPOPT maximum iterations
              
        solver.init()

        xMin = C.vertcat([ocp.variable(ocp.x[i].getName()).min.getValue() for i in range(ocp.x.size())])/stateScaling
        xMax = C.vertcat([ocp.variable(ocp.x[i].getName()).max.getValue() for i in range(ocp.x.size())])/stateScaling
     
        zMin = C.vertcat([ocp.variable(ocp.z[i].getName()).min.getValue() for i in range(ocp.z.size())])/algStateScaling
        zMax = C.vertcat([ocp.variable(ocp.z[i].getName()).max.getValue() for i in range(ocp.z.size())])/algStateScaling

        uMin = C.vertcat([ocp.variable(ocp.u[i].getName()).min.getValue() for i in range(ocp.u.size())])/controlScaling
        uMax = C.vertcat([ocp.variable(ocp.u[i].getName()).max.getValue() for i in range(ocp.u.size())])/controlScaling
        
        cMin = C.vertcat([ocp.variable(consList[i]).min.getValue() for i in range(len(consList))])/consScaling
        cMax = C.vertcat([ocp.variable(consList[i]).max.getValue() for i in range(len(consList))])/consScaling 
        
        solver.setInput(C.vertcat([NP.zeros(ocp.z.size()+ocp.x.size()),altUVals,cMin]),'lbg')	# g_L
        solver.setInput(C.vertcat([NP.zeros(ocp.z.size()+ocp.x.size()),altUVals,cMax]),'ubg')	# g_U        
  
  
        solver.setInput(C.vertcat([xMin,zMin,uMin]),  'lbx')					# u_L
        solver.setInput(C.vertcat([xMax,zMax,uMax]),  'ubx')					# u_U 
  

        solver.setInput(C.vertcat([x0,z0,u0]),'x0')   
  
        solver.evaluate()


        xzu = solver.output('x')
        split.setInput(xzu)
        split.evaluate()
        x0 = split.getOutput(0)
        z0 = split.getOutput(1)
        u0 = split.getOutput(2)        
        
        self.mSXF.setInput(x0,'x')
        self.mSXF.setInput(z0,'z')
        self.mSXF.setInput(u0,'p')
        
        self.mSXF.evaluate()
        
        y0 = self.mSXF.getOutput()        
        
        
        return x0*stateScaling,u0*controlScaling,z0*algStateScaling,y0*measurementScaling


    def computeYerr(self,openLoopControlTagsModelica,
                         openLoopControlTagsOlga,                         
                         measurementTagsModelica,
                         measurementTagsOlga,
                         measurementsMonitorModelica,
                         measurementsMonitorOlga,
                         controlTagsOlga,
                         dataFile ='yerrSim.mat',z0=None,u=None,x0=None):
        
        ocp = self.ocp
        if z0 is None:        
            z0 = C.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])
        if x0 is None:    
            x0 = C.vertcat([ocp.variable(ocp.x[k].getName()).initialGuess.getValue() for k in range(ocp.x.size())]) 
        if u is None:    
            u =  C.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])
                                                 
        data = sio.loadmat(dataFile)
        
        yh = data['yh']
        olgaTags = data['olgaTags']
        
        olgaTags = [str(ti).rstrip()  for ti in olgaTags]
        OlgaTagIndex = {}
        
        for k in range(len(olgaTags)):
            OlgaTagIndex[olgaTags[k]] = k
            
        
        OlgaF = {}
        for k in olgaTags:
            OlgaF[k] = interp1d([yi[OlgaTagIndex['TIME']] for yi in yh ],[yi[OlgaTagIndex[k]] for yi in yh ])
            
        
        finalTime = OlgaF['TIME'].x[-1]
        yFOlga = NP.array([OlgaF[ks](finalTime) for ks in measurementTagsOlga])
        yFMonitorOlga = NP.array([OlgaF[ks](finalTime) for ks in measurementsMonitorOlga])
        uFOlga = NP.array([OlgaF[ks](finalTime) for ks in controlTagsOlga])
        openUOlgaF = NP.array([OlgaF[ks](finalTime) for ks in openLoopControlTagsOlga])        
        
        
        
        xFM,uFM,zFM,yFM = self.findConstrainedSteadyState(openLoopControlTagsModelica,openUOlgaF,u,x0,z0)
        
        self.changeMeasurement(measurementsMonitorModelica);
        xFM,zFM,yFMonitorModelica = self.oneStep(xFM,uFM,zFM) 
        self.changeMeasurement(measurementTagsModelica);
        
        
        yerr = yFM - yFOlga 
        yerrMon = yFMonitorModelica - yFMonitorOlga
        uerr = uFM - uFOlga
        
        return yerr,yerrMon,uerr
        
    def initSimulation(self,x0):
        
        ocp = self.ocp

        ocp.variable('net.x0[1]').start = x0[7].getValue()
        ocp.variable('net.x0[2]').start = x0[8].getValue()
        ocp.variable('net.x0[3]').start = x0[9].getValue()
        ocp.variable('net.x0[4]').start = x0[10].getValue()
        ocp.variable('net.x0[5]').start = x0[11].getValue()
        ocp.variable('net.x0[6]').start = x0[12].getValue()
        ocp.variable('net.x0[7]').start = x0[13].getValue()
        ocp.variable('net.x0[8]').start = x0[14].getValue()
        ocp.variable('net.x0[9]').start = x0[15].getValue()
        ocp.variable('net.x0[10]').start = x0[16].getValue()

        w1_m_Ga = x0[7].getValue()
        w2_m_Ga = x0[10].getValue()
        
        R = 8.314472
        
        w1_T_a = ocp.variable('net.w1.par.T_a').start
        w1_M_G_a = ocp.variable('net.w1.par.M_G_a').start
        w1_V_a = ocp.variable('net.w1.par.V_a').start

        w2_T_a = ocp.variable('net.w2.par.T_a').start
        w2_M_G_a = ocp.variable('net.w2.par.M_G_a').start
        w2_V_a = ocp.variable('net.w2.par.V_a').start    
        
        PCA1Scaling = ocp.variable('net.w1.fin.p').nominal
        PCA2Scaling = ocp.variable('net.w2.fin.p').nominal
    
        ocp.variable('net.PCA1InitialSetPoint').start = R * w1_T_a * w1_m_Ga / (w1_M_G_a * w1_V_a) / PCA1Scaling
        
        ocp.variable('net.PCA2InitialSetPoint').start = R * w2_T_a * w2_m_Ga / (w2_M_G_a * w2_V_a) / PCA2Scaling
        
        self.ocp = ocp
        