# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 17:00:14 2014

@author: codas
"""

import casadi as C
import numpy as NP

class StateEstimator:
    
    
    ocp = None
    DT = None
    measurmentList = None
    mSX = None
    G = None 
    dxfdx0=None 
    pm = None
    pg = None
    Q = None
    R = None
        

    def initEstimator(self,ocp,DT,measuremntsList):
        
        


        sysIn = C.daeIn(x=ocp.x,z=ocp.z,p=ocp.u,t=ocp.t)
        sysOut = C.daeOut(ode=ocp.ode(ocp.x),alg=ocp.alg)
        odeF=C.SXFunction(sysIn,sysOut)
        odeF.addMonitor('inputs')
        odeF.addMonitor('outputs')        
        odeF.init()
        
        
        C.Integrator.loadPlugin("idas")
        G = C.Integrator("idas",odeF)					# Instantiate a CVodesIntegrator
#        f_d.setOption("abstol",self.CVODES_ABS_TOL)			# Absolute tolerance
#        f_d.setOption("reltol",self.CVODES_REL_TOL)			# Relative tolerance
#        f_d.setOption("max_multistep_order",15)  #for CVODES 
        G.setOption("max_multistep_order",5)  #for IDAS
        G.setOption("tf",DT)		        			# Integration length for each u^l
        
#==============================================================================
#        f_d.setOption('verbose',True)       
#        f_d.addMonitor('res')
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

        G.init()
        
        dxfdx0 = G.jacobian('x0','xf')
        dxfdx0.init()

        mSX = C.vertcat([self.ocp.variable(measuremntsList[k]).beq  for k in range(len(measuremntsList))])

        mSX.init()
        
        pm = mSX.jacobian()
        
        
        sysIn = C.daeIn(x=ocp.x,z=ocp.z,p=ocp.u,t=ocp.t)
        algF=C.SXFunction(sysIn,[ocp.alg])        
        algF.init()
        pg =algF.jacobian() 
        
        
        self.Q = NP.identity(ocp.x.size())
        self.R = NP.identity(len(measuremntsList))

        
        self.ocp = ocp
        self.DT = DT
        self.measurmentList = measuremntsList
        self.mSX = mSX
        self.G = G
        self.dxfdx0=dxfdx0 
        self.pm = pm
        self.pg = pg

    def DAE_Kalman(self,y,x0,u,z0):

        self.G.setInput(self.x0,'x0')        
        self.G.setInput(self.u,'p')
        self.G.setInput(self.z0,'z0')
        
        self.f_d.evaluate()
        
        x_bar = self.G.getOtput('xf')        
        z_bar = self.G.getOtput('zf')
        

        self.dxfdx0.setInput(self.X0, 'x0')
        self.dxfdx0.setInput(self.USOL[0], 'p')
        self.dxfdx0.setInput(self.Z0,'z0')
        self.dxfdx0.evaluate()
        dxfdx0 = self.dxfdx0.getOutput()        
        
        P_bar = C.mul(C.mul(dxfdx0,self.P_hat),dxfdx0.T) + self.Q
            
        ## estimated prediction
        
        self.mSX.setInput(x_bar, 'x')
        self.mSX.setInput(z_bar, 'z')
        self.mSX.setInput(u, 'p')
        self.mSX.evaluate()
        
        yP = self.mSX.getOutput()
        
        
        self.pm.setInput(x_bar,'x')
        self.pm.setInput(z_bar,'z')
        self.pm.setInput(u, 'p')
        self.pm.evaluate()
        
        self.pg.setInput(x_bar, 'x')
        self.pg.setInput(z_bar,'z')
        self.pg.setInput(u, 'p')
        self.pg.evaluate()
        
        HJx = self.pm.getOutput('xf') - self.pm.getOutput('zf') * NP.linalg.inv(self.pg.getOutput('zf')) * self.pg.getOutput('xf') 

        
        # Kalman Gain
        S = C.mul(C.mul(HJx,P_bar),HJx.T)+self.R
        SI = NP.linalg.inv(S)
        K = C.mul(C.mul(P_bar,HJx.T),SI)
    
        if NP.linalg.matrix_rank(NP.obsv(dxfdx0,HJx)) !=  self.ocp.x.size():
            AssertionError('Not observable')
    
        
        # Correction
        x_hat = x_bar + C.mul(K,y-yP)
        
        P_hat = C.mul((NP.eye(self.ocp.x.size())-C.mul(K,HJx)),P_bar)
        
        return {'x_hat':x_hat,'P_hat':P_hat,'yP':yP,'K':K}            
        
        
    def oneStep(self,x0,u,z0=None):
         
         
        outputSX = C.vertcat([self.ocp.variable(self.measurmentList[k]).beq  for k in range(len(self.measurmentList))])
        
        fout = C.SXFunction(C.daeIn(x=self.ocp.x,z=self.ocp.z,p=self.ocp.u,t=self.ocp.t),[outputSX])
        fout.init()
        
        if z0 is None:        
            z0 = C.vertcat([self.ocp.variable(self.ocp.z[k].getName()).start for k in range(self.ocp.z.size())])
        
        self.G.setInput(x0,'x0')
        self.G.setInput(z0,'z0')
        self.G.setInput(self.u,'p')
            
        self.G.evaluate()
            
        xk = self.G.getOutput('xf')
        zf = self.G.getOutput('zf')
                        
        return xk,zf           
        
        
        
        
        
        
        
        
        
        
        
        