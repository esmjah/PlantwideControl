# -*- coding: utf-8 -*-
"""
Created on Wed Jan 07 18:41:21 2015

@author: Administrator
"""
import scipy.io as sio
import numpy as NP

Ns = len(simTime)

nYkfMeas = len(measurementTagsModelica)
nYkfMoni = len(measurementsMonitorModelica)
nX = len(NP.array(xh[:][0]))
#nXs = len(NP.array(xSh[:][0]))

measurementTagsKf = NP.zeros((nYkfMeas,), dtype=NP.object)
measurementTagsKf[:] = measurementTagsModelica

measurementTagsMonitorKf = NP.zeros((nYkfMoni,), dtype=NP.object)
measurementTagsMonitorKf[:] = measurementsMonitorModelica

x = NP.zeros((Ns,nX))
#xs = NP.zeros((Ns,nXs))
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

yKfMeas = NP.zeros((Ns,nYkfMeas))
yKfMoni = NP.zeros((Ns,nYkfMoni))

yOlgaMeas = NP.zeros((Ns,nYolgaMeas))
yOlgaMoni = NP.zeros((Ns,nYolgaMoni))

u = NP.zeros((Ns,len(controlTagsModelica)))
Ju_NMPC = NP.zeros((Ns,len(controlTagsModelica)))
Ju1 = NP.zeros(Ns)
Ju2 = NP.zeros(Ns)
J = NP.zeros(Ns)
sim_time = NP.zeros(Ns)


for i in range(Ns):
    yKfMeas[:][i] = NP.array(yPh[:][i]).T
    yKfMoni[:][i] = NP.array(monitoredModelicah[:][i]).T
    yOlgaMeas[:][i] = NP.array(yh[:][i]).T
    yOlgaMoni[:][i] = NP.array(monitoredOlgah[:][i]).T
    x[:][i] = NP.array(xh[:][i]).T
    Ju_NMPC[:][i] = NP.array(Juh_NMPC[:][i]).T
    u[:][i] = NP.array(uh[:][i]).T
    J[i] = Objh_NMPC[i]
    sim_time[i] = simTime[i]
    
    
simData = {}
simData['yKfMeas'] = yKfMeas
simData['yKfMoni'] = yKfMoni
simData['yOlgaMeas'] = yOlgaMeas
simData['yOlgaMoni'] = yOlgaMoni
simData['x'] = x
simData['u_NMPC'] = u
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
simData['J_NMPC'] = J
simData['Ju_NMPC'] = Ju_NMPC
simData['simTime'] = sim_time

sio.savemat('SavedResults\\D1D2_NMPC.mat', simData)
