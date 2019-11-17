# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 17:00:14 2014

@author: codas
"""

import casadi as C
import numpy as NP
import scipy


def obsv(A, C):       


    # Convert input parameters to matrices (if they aren't already)
    amat = NP.mat(A)
    cmat = NP.mat(C)
    n = NP.shape(amat)[0]

    # Construct the controllability matrix
    obsv = cmat
    for i in range(1, n):
        obsv = NP.vstack((obsv, cmat*amat**i))
         
    return obsv

    
def kalmd(A,C,Q,R):
 
 
    #first, try to solve the ricatti equation
    X = NP.matrix(scipy.linalg.solve_discrete_are(A.T, C.T, Q, R))
     
    #compute the LQR gain
     
    L = (A*X*C.T) * scipy.linalg.inv(C*X*C.T+R)
     
    eigVals, eigVecs = scipy.linalg.eig(A-L*C)
     
    return X,eigVals,L     
    
class StateEstimator:
    
    
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

    def initEstimator(self,ocp,DT,measuremntsList,Q=None,R=None):
        
           
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
        G.setOption("reltol",1e-6)  #for IDAS
        G.setOption("abstol",1e-6)  #for IDAS 
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

        if Q is None:
            self.Q = NP.identity(ocp.x.size())
        else:
            self.Q = Q            
            
        if R is None:            
            self.R = NP.identity(len(measuremntsList))
        else:
            self.R = R
            
        self.P_hat = NP.identity(ocp.x.size())

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
        
        
    def kalmanLinear(self,y,u,x0=None,z0=None):
        
        if self.algStateScaling.size() == 0:
            out = self.ODE_kalmanLinear(y,u,x0)
        else:
            error('not Implmented')
        
        return out
        
        
    def computeKalmanSystem(self,u,x0):

        if x0 is None:
            x0 = self.x0
        else:
            x0 = NP.copy(x0)/self.stateScaling
                                    
        u = u/self.controlScaling                                    
                                    
        self.dxfdx0.setInput(x0, 'x0')
        self.dxfdx0.setInput(u, 'p')
        self.dxfdx0.evaluate()
        A = NP.array(self.dxfdx0.getOutput(0))

        
        self.pm.setInput(x0,'x')
        self.pm.setInput(u, 'p')
        self.pm.evaluate()
        C = NP.array(self.pm.getOutput(0))     
        D = self.pm.getOutput(3)      
        y0 = self.pm.getOutput(2)
        
        X,eigVals,L = kalmd(A,C,self.Q,self.R)          
        
        self.dxfdu.setInput(x0, 'x0')
        self.dxfdu.setInput(u, 'p')
        self.dxfdu.evaluate()
        B = self.dxfdu.getOutput(0)      
        
        self.AK = A-L*NP.array(C)
        self.BK = B-L*NP.array(D)       
        self.CK = C  
        self.DK = D  
        self.L = L
        self.eigVals = eigVals
        self.x_ref = x0
        self.u_ref = u
        self.y_ref = y0
        
        
        
    def ODE_kalmanLinear(self,y,u,xk_1):
        
        u = u/self.controlScaling        
        
        dxk_1 = NP.copy(xk_1)/self.stateScaling - self.x_ref
        dy = NP.copy(y)/self.measurementScaling - self.y_ref
        du = u - self.u_ref
        
        dx_k = C.mul(self.AK , dxk_1) + C.mul(self.BK , du) + C.mul(self.L , dy)
        
        x_k = self.x_ref + dx_k
       
        y_k = self.y_ref + C.mul(self.CK , dx_k)  
       
       
        return {'x_hat':x_k*self.stateScaling,'yP':y_k*self.measurementScaling,'zh':self.algStateScaling}            
        
       
        

    def kalman(self,y,u,x0=None,z0=None):
        
        if self.algStateScaling.size() == 0:
            out = self.ODE_kalman(y,u,x0)
        else:
            out = self.DAE_Kalman(y,u,x0,z0)
        
        return out


    def ODE_kalman(self,y,u,x0=None):


        if x0 is None:    
            x0 = self.x0
        else:
            x0 = NP.copy(x0)/self.stateScaling
            
        y = NP.copy(y)/self.measurementScaling      
        u = NP.copy(u)/self.controlScaling
        
        self.dxfdx0.setInput(x0, 'x0')
        self.dxfdx0.setInput(u, 'p')
        self.dxfdx0.evaluate()
        dxfdx0 = self.dxfdx0.getOutput(0)
        
        x_bar = self.dxfdx0.getOutput('xf')
        
        P_bar = C.mul(C.mul(dxfdx0,self.P_hat),dxfdx0.T) + self.Q
            
        ## estimated prediction            
        
        self.pm.setInput(x_bar,'x')
        self.pm.setInput(u, 'p')
        self.pm.evaluate()

        yP = self.pm.getOutput(2)        
           
        self.pg.setInput(x_bar, 'x')
        self.pg.setInput(u, 'p')
        self.pg.evaluate()
        
         
        
        HJx = self.pm.getOutput(0) 

#        print HJx,P_bar,self.R
        # Kalman Gain
        S = C.mul(C.mul(HJx,P_bar),HJx.T)+self.R
        SI = NP.linalg.inv(S)
        K = C.mul(C.mul(P_bar,HJx.T),SI)
    
        if NP.linalg.matrix_rank(obsv(dxfdx0,HJx)) !=  self.ocp.x.size():
            AssertionError('Not observable')
    
        
        # Correction
        x_c =  C.mul(K,y-yP) 
        
        
        x_hat = x_bar + x_c
        #z_hat = self.findZ(u,x0=x_hat*self.stateScaling,z0=z_bar*self.algStateScaling)/self.algStateScaling
        
        self.mSXF.setInput(x_hat,'x')        
        self.mSXF.setInput(u,'p')
        
        self.mSXF.evaluate()

        yC = self.mSXF.getOutput()  

#        print 'Kalman Correction',x_c                
#        print 'scaled   Measurment',y        
#        print 'predictedMeasurment',yP
#        print 'correctedMeasurment',yC        
        
        P_hat = C.mul((NP.eye(self.ocp.x.size())-C.mul(K,HJx)),P_bar)
        self.P_hat = P_hat;
        
        self.x0 = NP.copy(x_hat)
        
        return {'x_hat':x_hat*self.stateScaling,'P_hat':P_hat,'yC':yC*self.measurementScaling,'yP':yP*self.measurementScaling,'K':K,'zh':self.algStateScaling}            
        

    def DAE_Kalman(self,y,u,x0=None,z0=None):
        
        if z0 is None:        
            self.findZ(u)
            z0 = self.z0
        else:
            z0 = NP.copy(z0)/self.algStateScaling
                
        if x0 is None:    
            x0 = self.x0
        else:
            x0 = NP.copy(x0)/self.stateScaling
            
        y = NP.copy(y)/self.measurementScaling      
        u = NP.copy(u)/self.controlScaling
        
        self.dxfdx0.setInput(x0, 'x0')
        self.dxfdx0.setInput(u, 'p')
        self.dxfdx0.setInput(z0,'z0')
        self.dxfdx0.evaluate()
        dxfdx0 = self.dxfdx0.getOutput(0)
        
        x_bar = self.dxfdx0.getOutput('xf')
        z_bar = self.dxfdx0.getOutput('zf')          
        
        P_bar = C.mul(C.mul(dxfdx0,self.P_hat),dxfdx0.T) + self.Q
            
        ## estimated prediction
              
        
        self.pm.setInput(x_bar,'x')
        self.pm.setInput(z_bar,'z')
        self.pm.setInput(u, 'p')
        self.pm.evaluate()

        yP = self.pm.getOutput(2)        
        

        
        
        self.pg.setInput(x_bar, 'x')
        self.pg.setInput(z_bar,'z')
        self.pg.setInput(u, 'p')
        self.pg.evaluate()
        
         
        
        HJx = self.pm.getOutput(0) - C.mul(C.mul(self.pm.getOutput(1), NP.linalg.inv(self.pg.getOutput(1))) , self.pg.getOutput(0)) 

#        print HJx,P_bar,self.R
        # Kalman Gain
        S = C.mul(C.mul(HJx,P_bar),HJx.T)+self.R
        SI = NP.linalg.inv(S)
        K = C.mul(C.mul(P_bar,HJx.T),SI)
    
        if NP.linalg.matrix_rank(obsv(dxfdx0,HJx)) !=  self.ocp.x.size():
            AssertionError('Not observable')
    
        
        # Correction
        x_c =  C.mul(K,y-yP) 
        
        
        x_hat = x_bar + x_c
        z_hat = self.findZ(u*self.controlScaling,x0=x_hat*self.stateScaling,z0=z_bar*self.algStateScaling)/self.algStateScaling
        z_hat = self.algStateScaling
        
        self.mSXF.setInput(x_hat,'x')        
        self.mSXF.setInput(z_hat,'z')
        self.mSXF.setInput(u,'p')
        
        self.mSXF.evaluate()

        yC = self.mSXF.getOutput()  

#        print 'Kalman Correction',x_c                
#        print 'scaled   Measurment',y        
#        print 'predictedMeasurment',yP
#        print 'correctedMeasurment',yC        
        
        P_hat = C.mul((NP.eye(self.ocp.x.size())-C.mul(K,HJx)),P_bar)
        self.P_hast = P_hat
        
        self.z0 = NP.copy(z_hat)  ## this is just an estimation
        self.x0 = NP.copy(x_hat)
        
        return {'x_hat':x_hat*self.stateScaling,'P_hat':P_hat,'yC':yC*self.measurementScaling,'yP':yP*self.measurementScaling,'K':K,'zh':self.z0*self.algStateScaling}            
        
    def getState(self):
        return self.x0*self.stateScaling
        
    def setState(self,x0,z0):
        self.x0 = x0/self.stateScaling
        self.z0 = z0/self.algStateScaling
        
        
    def oneStep(self,x0,u,z0=None):
        
        if z0 is None:        
            z0 = self.z0
        else:
            z0 = z0/self.algStateScaling
            
        if x0 is None:    
            x0 = self.x0
        else:
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