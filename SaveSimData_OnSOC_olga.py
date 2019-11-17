# -*- coding: utf-8 -*-
"""
Created on Wed Jan 07 18:41:21 2015

@author: Administrator
"""
import scipy.io as sio
import numpy as NP

NIT = len(uh)

nYkfMeas = len(measurementTagsModelica)
nYkfMoni = len(measurementsMonitorModelica)
nX = len(NP.array(xh[:][0]))
# nXs = len(NP.array(xSh[:][0]))

measurementTagsKf = NP.zeros((nYkfMeas,), dtype=NP.object)
measurementTagsKf[:] = measurementTagsModelica

measurementTagsMonitorKf = NP.zeros((nYkfMoni,), dtype=NP.object)
measurementTagsMonitorKf[:] = measurementsMonitorModelica

x = NP.zeros((NIT, nX))
# xs = NP.zeros((NIT,nXs))
states = NP.zeros((nX,), dtype=NP.object)
xMin = NP.zeros(nX)
xMax = NP.zeros(nX)
xOptimal = NP.zeros(nX)
for k in range(nX):
    states[k] = ocp.x[k].getName()
    xMin[k] = ocp.variable(ocp.x[k].getName()).min.getValue()
    xMax[k] = ocp.variable(ocp.x[k].getName()).max.getValue()
    xOptimal[k] = NP.array(xOpt[k])[0][0]

nYolgaMeas = len(measurementTagsOlga)
nYolgaMoni = len(measurementsMonitorOlga)

measurementTags = NP.zeros((nYolgaMeas,), dtype=NP.object)
measurementTags[:] = measurementTagsOlga

measurementTagsMonitor = NP.zeros((nYolgaMoni,), dtype=NP.object)
measurementTagsMonitor[:] = measurementsMonitorOlga

nU = len(controlTagsModelica)
uMin = NP.zeros(nU)
uMax = NP.zeros(nU)
uOptimal = NP.zeros(nU)
for k in range(nU):
    uMin[k] = ocp.variable(ocp.u[k].getName()).min.getValue()
    uMax[k] = ocp.variable(ocp.u[k].getName()).max.getValue()
    uOptimal[k] = NP.array(uOpt[k])[0][0]

yOptimal = NP.zeros(len(measurementTagsModelica))
for k in range(len(measurementTagsModelica)):
    yOptimal[k] = NP.array(yOpt[k])[0][0]

measurementTagsOptm = NP.zeros((len(measurementTagsModelica),), dtype=NP.object)
measurementTagsOptm[:] = measurementTagsModelica

controlTagsKf = NP.zeros((nU,), dtype=NP.object)
controlTagsKf[:] = controlTagsModelica

controlTags = NP.zeros((len(controlTagsOlga[0:nU]),), dtype=NP.object)
controlTags[:] = controlTagsOlga[0:nU]

yKfMeas = NP.zeros((NIT, nYkfMeas))
yKfMoni = NP.zeros((NIT, nYkfMoni))

yOlgaMeas = NP.zeros((NIT, nYolgaMeas))
yOlgaMoni = NP.zeros((NIT, nYolgaMoni))

u = NP.zeros((NIT, len(controlTagsModelica)))
Ju_OnSOC = NP.zeros((NIT, len(controlTagsModelica)))
Ju1 = NP.zeros(NIT)
Ju2 = NP.zeros(NIT)
J = NP.zeros(NIT)

for i in range(NIT):
    yKfMeas[:][i] = NP.array(yPh[:][i]).T
    yKfMoni[:][i] = NP.array(monitoredModelicah[:][i]).T
    yOlgaMeas[:][i] = NP.array(yh[:][i]).T
    yOlgaMoni[:][i] = NP.array(monitoredOlgah[:][i]).T
    x[:][i] = NP.array(xh[:][i]).T
    Ju_OnSOC[:][i] = NP.array(Juh_OnSOC[:][i]).T
    u[:][i] = NP.array(uh[:][i]).T
    J[i] = Objh[i]

simData = {}
simData['yKfMeas'] = yKfMeas
simData['yKfMoni'] = yKfMoni
simData['yOlgaMeas'] = yOlgaMeas
simData['yOlgaMoni'] = yOlgaMoni
simData['x'] = x
simData['u'] = u
simData['measurementTagsKf'] = measurementTagsKf
simData['measurementTagsKf'] = measurementTagsKf
simData['measurementTagsOpt'] = measurementTagsOptm
simData['measurementTagsMonitorKf'] = measurementTagsMonitorKf
simData['measurementTags'] = measurementTags
simData['measurementTagsMonitor'] = measurementTagsMonitor
simData['controlTagsKf'] = controlTagsKf
simData['controlTags'] = controlTags
simData['states'] = states
simData['uMin'] = uMin
simData['uMax'] = uMax
simData['xMin'] = xMin
simData['xMax'] = xMax
simData['uOptimal'] = uOptimal
simData['xOptimal'] = xOptimal
simData['yOptimal'] = yOptimal
simData['J'] = J
simData['Ju_OcSOC'] = Ju_OnSOC

sio.savemat('SimData_OnSOC_Nominal_2018.06.25_Olga.mat', simData)
