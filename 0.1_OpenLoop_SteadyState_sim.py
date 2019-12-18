# -*- coding: utf-8 -*-
# Authors: Andres Codas, Esmaeil Jahanshahi
# Norwegian University of Science and Technology

# This script calculates the sates at the stationary point (x*) of the open-loop system for given u, where f(x*,u)=0
# The x* is used to initialize the dynamic simulations
# Also, we get the liner state-space model (A, B, C, D) of the open-loop
# The liner state-space model (A, B, C, D) is used for controllability analysis in the regulatory control layer
# This scripts uses "ocpGaslift_OpenLoop_Sim.mop" as the model. This model takes the 4 disturbances as inputs
# The model has 9 inputs: [Alpha1, Alpha2, P_res1, P_res2, inputGas1, inputGas2, zManifold1, zManifold2, zPipe]

import numpy as np
import os
import sys

_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(_path)
sys.path.append(_path + '\\modelCompiler')
sys.path.append(_path + '\\DaeModel')

import modelCompiler
import DaeModel

reload(modelCompiler)
reload(DaeModel)

measurementTagsModelica = ['net.p.w_mix_out',
                           'net.w1sim.w_L_out',
                           'net.w2sim.w_L_out',
                           'net.p.fin.p',
                           'net.p.Alpha_L1_av',
                           'net.w1sim.P_bh',
                           'net.w2sim.P_bh',
                           'net.w1sim.fin.p','net.w2sim.fin.p','net.p.deltaP_valve','net.w1sim.deltaP_valve','net.w2sim.deltaP_valve']

models_list = modelCompiler.getModelList('.', True)
ocp = modelCompiler.getOCP('ocpGaslift_OpenLoop_Sim', models_list)

ocp.makeSemiExplicit()
ocp.eliminateIndependentParameters()
ocp.eliminateDependentParameters()
ocp.eliminateDependentParameterInterdependencies()
ocp.eliminateAlgebraic()

# u = C.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])

u = np.array([0, 0, 160, 170, 1, 1, 0.5, 0.5, 0.5])  # Initial conditions
# u = np.array([0, 0, 160, 170, 1.29622, 1.32441, 0.580464, 0.608216, 0.576931])  # Steady-state optimal

daeModel = DaeModel.DaeModel()
DT = 10
daeModel.init(ocp, DT, measurementTagsModelica)

x0, z0, y0 = daeModel.findSteadyState(u,None,None,0,consList=['net.w1sim.z1','net.w2sim.z1','net.p.z'])
measurementTags = ['net.w1sim.fin.p', 'net.Dp_w1', 'net.w2sim.fin.p', 'net.Dp_w2', 'net.p.fin.p', 'net.p.deltaP_valve']

a, b, c, d = daeModel.computeJacobians(x0, u, measurementTags)

#sio.savemat('C:\Users\Admin\Dropbox\Journal_Papers\Plantwide\SensitivityAnalysis\jacobians_struc2_zOpt', mdict={'A': np.matrix(a),'B': np.matrix(b),'C': np.matrix(c),'D': np.matrix(d)})

print 'y0', y0
print 'x0', x0

