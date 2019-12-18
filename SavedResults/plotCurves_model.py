import matplotlib.pyplot as plt
plt.close('all')

setPointsModelicaTags = ['uTOP', 'uWHD1', 'uWHD2', 'net.PCA1.setPoint', 'net.PCA2.setPoint']
controledVariableModelicaTags = ['net.p.deltaP_valve', 'net.w1.deltaP_valve', 'net.w2.deltaP_valve',
                                 'net.PCA1.measurement', 'net.PCA2.measurement']

setPointIndexInControl = []
openLoop = False

k0 = 0
tc = range(0, len(uh))
tc = [(k0 + ti) * DT for ti in tc]

nInputs = uh[0].size()

for si in range(nInputs):
   ui = [i[si] for i in uh]
   plt.figure()
   plt.plot(tc,ui,label = ocp.u[si].getName())
   plt.legend()

for tag in setPointsModelicaTags:
    for si in range(nInputs):
        if ocp.u[si].getName() == tag:
            setPointIndexInControl.append(si)

# plot outputs and predicted outputs
t = range(1, len(yPh) + 1)
t = [(k0 + ti) * DT for ti in t]

nSignals = len(measurementTagsSim)
for si in range(nSignals):

    yiM = [i[si] for i in yPh]
    yiN = [yh[i][si] for i in range(len(yh))]

    plt.figure()
    ploted = False
    if not openLoop:
        if controlTagsSim:
            for cvI in range(len(controledVariableModelicaTags)):
                if controledVariableModelicaTags[cvI] == measurementTagsModelica[si]:
                    ui = [i[setPointIndexInControl[cvI]] for i in uh]

                    plt.plot(t, yiN, t, yiM, tc, ui)
                    plt.legend([measurementTagsSim[si], measurementTagsModelica[si],
                                ocp.u[setPointIndexInControl[cvI]].getName()])
                    ploted = True

            if not ploted:
                plt.plot(t, yiN, t, yiM, 'r')
                plt.legend([measurementTagsSim[si], measurementTagsModelica[si]])
                ploted = True

    if not ploted:
        for cvI in range(len(controledVariableModelicaTags)):
            if controledVariableModelicaTags[cvI] == measurementTagsModelica[si]:
                yiO = [SimF[measurementTagsSim[si]](t[i]) for i in range(len(yPh))]
                ui = [i[setPointIndexInControl[cvI]] for i in uh]

                plt.plot(t, yiN, t, yiO, t, yiM, tc, ui)
                plt.legend([measurementTagsSim[si], measurementTagsSim[si], measurementTagsModelica[si],
                            ocp.u[setPointIndexInControl[cvI]].getName()])
                ploted = True

    if not ploted:
        yiO = [SimF[measurementTagsSim[si]](t[i]) for i in range(len(yPh))]
        plt.plot(t, yiN, t, yiO, t, yiM)
        plt.legend([measurementTagsSim[si], measurementTagsSim[si], measurementTagsModelica[si]])
        ploted = True

        # plt.figure()
    # yerr = [SimF[measurementTagsSim[si]](t[i])-yPh[i][si] for i in range(len(yPh)) ]
    # plt.plot(t,yerr, label = 'Error ' + measurementTagsSim[si])
    plt.legend()
k0 = 0
t = range(len(xh))
t = [(k0 + ti) * DT for ti in t]

nStates = xh[0].size()
xSh = []
plotStates = True
if plotStates:
    for si in range(nStates):
        xi = [i[si] for i in xh]
        plt.figure()
        plt.plot(t, xi, label=ocp.x[si].getName())
        plt.legend()
        if len(xSh) > 0:
            xih = [i[si] for i in xSh]
            plt.plot(t, xih[:len(t)])
            plt.legend()

        plt.plot(t, [xOpt[si] for ti in t])
        plt.plot(t, [xMin[si] for ti in t], "--rx")
        plt.plot(t, [xMax[si] for ti in t], "--rx")


if len(monitoredSimh) > 0:
    t = range(len(monitoredSimh))
    t = [(k0 + ti + 1) * DT for ti in t]
    for si in range(len(measurementsMonitorModelica)):
        plt.figure()
        v = [monitoredSimh[k][si] for k in range(len(monitoredSimh))]
        v2 = [monitoredModelicah[k][si] for k in range(len(monitoredModelicah))]
        plt.plot(t, v, t, v2)
        plt.legend([measurementsMonitorSim[si], measurementsMonitorModelica[si]])
