# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 17:00:14 2014

@author: codas
"""

import casadi as C
import numpy as NP
import scipy


     
    
class ModelSimulator:
    
    
    ocp = None
    DT = None
    measurmentList = None
    mSXF = None
    G = None 
    dxfdx0=None 
    pm = None
    pg = None
    Q = None
    R = None
    P_hat = None    
    
    x0 = None
    z0 = None
    measurementScaling = None
    stateScaling = None
    algStateScaling = None
    controlScaling = None
    
    algF = None # remove!

    def initSimulator(self,ocp,DT,measuremntsList):
        
           
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
        G.setOption("max_step_size",1)
        G.setOption("tf",DT)		        			
        G.init()
#==============================================================================
#        G.setOption('verbose',True)       
#        G.addMonitor('res')
#        f_d.addMonitor('djacB')
#        f_d.addMonitor('inputs')        
#        f_d.addMonitor('outputs')        
#        f_d.addMonitor('bjacB')
#        f_d.addMonitor('jtimesB')
#        f_d.addMonitor('psetup')
#        f_d.addMonitor('psetupB')
#        f_d.addMonitor('psolveB')
#        f_d.addMonitor('resB')
#        f_d.addMonitor('resS')
#        f_d.addMonitor('rhsQB')
#==============================================================================


        
        dxfdx0 = G.jacobian('x0','xf')
        dxfdx0.init()

        dxfdu = G.jacobian('p','xf')
        dxfdu.init()

        
        mSX = C.vertcat([ocp.variable(measuremntsList[k]).beq  for k in range(len(measuremntsList))])
        mSX = C.substitute(mSX,C.vertcat([ocp.x,ocp.z,ocp.u]),C.vertcat([stateScaling*ocp.x,algStateScaling*ocp.z,controlScaling*ocp.u]))/measurementScaling
        mSXF = C.SXFunction(sysIn,[mSX])
        mSXF.init()
        
        pm = mSXF.jacobian()
        pm.init()
        
        d_alg_d_x = odeF.jac('x','alg')
        d_alg_d_z = odeF.jac('z','alg')
        pg = C.SXFunction(sysIn,[d_alg_d_x, d_alg_d_z,algS])
        pg.init()
        
        d_m_d_x = mSXF.jac('x',0)
        d_m_d_z = mSXF.jac('z',0)
        d_m_d_u = mSXF.jac('p',0)
        
        pm = C.SXFunction(sysIn,[d_m_d_x, d_m_d_z,mSX,d_m_d_u])            
        pm.init()
                
        
        self.z0 = C.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])/algStateScaling
        self.x0 = C.vertcat([ocp.variable(ocp.x[k].getName()).initialGuess.getValue() for k in range(ocp.x.size())])/stateScaling

            
        self.measurementScaling = measurementScaling
        self.stateScaling = stateScaling
        self.algStateScaling = algStateScaling
        self.controlScaling = controlScaling

        self.ocp = ocp
        self.DT = DT
        self.measurmentList = measuremntsList
        self.mSXF = mSXF
        self.G = G
        self.dxfdx0=dxfdx0
        self.dxfdu =dxfdu
        self.pm = pm
        self.pg = pg
        
        
    def Model_sim(self,u,x0=None):

        if x0 is None:    
            x0 = self.x0
        else:
            x0 = NP.copy(x0)/self.stateScaling
               
        u = NP.copy(u)/self.controlScaling
        
        self.dxfdx0.setInput(x0, 'x0')
        self.dxfdx0.setInput(u, 'p')
        self.dxfdx0.evaluate()
        
        x_k = self.dxfdx0.getOutput('xf')

        ## estimated prediction            
        
        self.pm.setInput(x_k,'x')
        self.pm.setInput(u, 'p')
        self.pm.evaluate()

        y_k = self.pm.getOutput(2)        
           
        self.x0 = NP.copy(x_k)
        
        return {'x_k':x_k*self.stateScaling,'y_k':y_k*self.measurementScaling}            
        
    def setState(self,x0,z0):
        self.x0 = x0/self.stateScaling
        self.z0 = z0/self.algStateScaling
        
    def findSteadyState(self,u,x0=None,z0=None,simCount=0,consList=[]):
        
        if z0 is not None:        
            z0 = z0/self.algStateScaling
        else:
            z0 = self.z0
        if x0 is not None:    
            x0 = x0/self.stateScaling
        else:
            x0 = self.x0
                
        ocp = self.ocp

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
        solver.setOption('print_user_options','yes')
        solver.setOption("max_iter",200)	# IPOPT maximum iterations
              
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