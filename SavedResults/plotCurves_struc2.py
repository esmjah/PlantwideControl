import matplotlib.pyplot as plt
plt.close('all')

setPointsModelicaTags = ['uTOP', 'uWHD1', 'uWHD2', 'net.PCA1.setPoint', 'net.PCA2.setPoint']
controledVariableModelicaTags = ['net.p.deltaP_valve', 'net.w1.deltaP_valve', 'net.w2.deltaP_valve',
                                 'net.PCA1.measurement', 'net.PCA2.measurement']

setPointIndexInControl = []

nInputs = uh[0].size()
for tag in setPointsModelicaTags:
    for si in range(nInputs):
        if ocp.u[si].getName() == tag:
            setPointIndexInControl.append(si)

k0 = 0
t = range(len(xh))
t = [(k0 + ti) * DT for ti in t]

nStates = xh[0].size()

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

tc = range(0, len(uh))
tc = [(k0 + ti) * DT for ti in tc]

# nInputs = uh[0].size()
#
# for si in range(nInputs):
#
#    ui = [i[si] for i in uh]
#    plt.figure()
#    plt.plot(tc,ui,label = ocp.u[si].getName())
#    plt.legend()


# plot outputs and predicted outputs
t = range(1, len(yPh) + 1)
t = [(k0 + ti) * DT for ti in t]

nSignals = len(measurementTagsOlga)
for si in range(nSignals):

    yiM = [i[si] for i in yPh]
    yiN = [yh[i][si] for i in range(len(yh))]

    plt.figure()
    ploted = False
    if not openLoop:
        if controlOlga:
            for cvI in range(len(controledVariableModelicaTags)):
                if controledVariableModelicaTags[cvI] == measurementTagsModelica[si]:
                    ui = [i[setPointIndexInControl[cvI]] for i in uh]

                    plt.plot(t, yiN, t, yiM, tc, ui)
                    plt.legend([measurementTagsOlga[si], measurementTagsModelica[si],
                                ocp.u[setPointIndexInControl[cvI]].getName()])
                    ploted = True

            if not ploted:
                plt.plot(t, yiN, t, yiM, 'r')
                plt.legend([measurementTagsOlga[si], measurementTagsModelica[si]])
                ploted = True

    if not ploted:
        for cvI in range(len(controledVariableModelicaTags)):
            if controledVariableModelicaTags[cvI] == measurementTagsModelica[si]:
                yiO = [OlgaF[measurementTagsOlga[si]](t[i]) for i in range(len(yPh))]
                ui = [i[setPointIndexInControl[cvI]] for i in uh]

                plt.plot(t, yiN, t, yiO, t, yiM, tc, ui)
                plt.legend([measurementTagsOlga[si], measurementTagsOlga[si], measurementTagsModelica[si],
                            ocp.u[setPointIndexInControl[cvI]].getName()])
                ploted = True

    if not ploted:
        yiO = [OlgaF[measurementTagsOlga[si]](t[i]) for i in range(len(yPh))]
        plt.plot(t, yiN, t, yiO, t, yiM)
        plt.legend([measurementTagsOlga[si], measurementTagsOlga[si], measurementTagsModelica[si]])
        ploted = True

        # plt.figure()
    # yerr = [OlgaF[measurementTagsOlga[si]](t[i])-yPh[i][si] for i in range(len(yPh)) ]
    # plt.plot(t,yerr, label = 'Error ' + measurementTagsOlga[si])
    plt.legend()

if len(monitoredOlgah) > 0:
    t = range(len(monitoredOlgah))
    t = [(k0 + ti + 1) * DT for ti in t]
    for si in range(len(measurementsMonitorModelica)):
        plt.figure()
        v = [monitoredOlgah[k][si] for k in range(len(monitoredOlgah))]
        v2 = [monitoredModelicah[k][si] for k in range(len(monitoredModelicah))]
        plt.plot(t, v, t, v2)
        plt.legend([measurementsMonitorOlga[si], measurementsMonitorModelica[si]])
