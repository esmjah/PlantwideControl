import casadi as C
import numpy as NP
import matplotlib.pyplot as plt
import ConfigParser
import sys


class SingleShooting:
    ocp = None
    ###### Optimizer Tolerances
    INTG_ABS_TOL = 1e-6  # Integrator absolute tolerance
    INTG_REL_TOL = 1e-6  # Integrator relative tolerance
    IDAS_MAX_STEP_SIZE = 2  # IDAS max step size

    IPOPT_tol = 1e-3  ## master tolerance 1e-3
    IPOPT_dual_inf_tol = 1e-2  ## dual infeasibility
    IPOPT_constr_viol_tol = 1e-3  ## primal infeasibility 1e-3
    IPOPT_compl_inf_tol = 1e-2  ## complementarity
    IPOPT_print_level = 5
    IPOPT_max_iter = 30
    IPOPT_max_cpu_time = 60 * 3600 * 24  ## one day of computation should be enough right? :)

    n_x = None
    n_u = None
    n_v = None
    n_k = None

    DT = None

    V = None
    U = None

    X0 = None

    USOL = None
    solver = None
    nlp = None

    odeF = None
    G = None
    controlCost = 0

    xOpt = None
    uOpt = None
    plotTags = []
    plotDic = {}

    def defineOCP(self, ocp, DT=20, controlCost=0, xOpt=[], uOpt=[], finalStateCost=0.0, deltaUCons=[]):

        self.ocp = ocp
        ocp = self.ocp
        self.DT = DT
        self.n_k = int(self.ocp.tf / self.DT)
        self.controlCost = controlCost

        stateScaling = C.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
        algStateScaling = C.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
        controlScaling = C.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])

        xOpt = xOpt / stateScaling
        uOpt = uOpt / controlScaling
        self.xOpt = xOpt
        self.uOpt = uOpt

        self.stateScaling = C.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
        self.algStateScaling = C.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
        self.controlScaling = C.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])

        odeS = C.substitute(ocp.ode(ocp.x), C.vertcat([ocp.x, ocp.z, ocp.u]), C.vertcat(
            [stateScaling * ocp.x, algStateScaling * ocp.z, controlScaling * ocp.u])) / stateScaling
        algS = C.substitute(ocp.alg, C.vertcat([ocp.x, ocp.z, ocp.u]),
                            C.vertcat([stateScaling * ocp.x, algStateScaling * ocp.z, controlScaling * ocp.u]))
        ltermS = C.substitute(ocp.lterm, C.vertcat([ocp.x, ocp.z, ocp.u]),
                              C.vertcat([stateScaling * ocp.x, algStateScaling * ocp.z, controlScaling * ocp.u]))

        sysIn = C.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t)
        sysOut = C.daeOut(ode=odeS, alg=algS, quad=ltermS)
        odeF = C.SXFunction(sysIn, sysOut)
        odeF.init()

        C.Integrator.loadPlugin("idas")
        G = C.Integrator("idas", odeF)
        G.setOption("reltol", self.INTG_REL_TOL)  # for CVODES and IDAS
        G.setOption("abstol", self.INTG_ABS_TOL)  # for CVODES and IDAS
        G.setOption("max_multistep_order", 5)  # for CVODES and IDAS
        G.setOption("max_step_size", self.IDAS_MAX_STEP_SIZE)  # for IDAS only
        G.setOption("tf", self.DT)
        self.G = G

        # ==============================================================================
        #        G.setOption('verbose',True)
        #        G.addMonitor('res')
        #        G.addMonitor('inputs')
        #        G.addMonitor('outputs')
        # G.addMonitor('djacB')
        #         G.addMonitor('bjacB')
        #         G.addMonitor('jtimesB')
        #         G.addMonitor('psetup')
        #         G.addMonitor('psetupB')
        #         G.addMonitor('psolveB')
        #         G.addMonitor('resB')
        #         G.addMonitor('resS')
        #         G.addMonitor('rhsQB')
        # ==============================================================================
        G.init()

        self.n_u = self.ocp.u.size()
        self.n_x = self.ocp.x.size()
        self.n_v = self.n_u * self.n_k

        self.V = C.MX.sym("V", int(self.n_v), 1)
        self.U = self.splitVariables(self.V)

        uMin = C.vertcat(
            [self.ocp.variable(self.ocp.u[i].getName()).min.getValue() for i in range(self.n_u)]) / controlScaling
        uMax = C.vertcat(
            [self.ocp.variable(self.ocp.u[i].getName()).max.getValue() for i in range(self.n_u)]) / controlScaling
        UMIN = C.vertcat([uMin for k in range(self.n_k)])
        UMAX = C.vertcat([uMax for k in range(self.n_k)])

        if len(deltaUCons) > 0:
            addDeltaUCons = True
            deltaUCons = deltaUCons / self.controlScaling
        else:
            addDeltaUCons = False
        ## PATH CONSTRAINTS

        pathIn = C.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t)

        pathVarNames = [sv.getName() for sv in ocp.beq(ocp.path)]
        pathScaling = C.vertcat([ocp.nominal(pv) for pv in pathVarNames])

        pathS = C.substitute(ocp.beq(ocp.beq(ocp.path)), C.vertcat([ocp.x, ocp.z, ocp.u]), C.vertcat(
            [stateScaling * ocp.x, algStateScaling * ocp.z, controlScaling * ocp.u])) / pathScaling
        pathConstraints = C.SXFunction(pathIn, [pathS]);

        pathMax = C.vertcat(
            [ocp.variable(pathVarNames[i]).max.getValue() for i in range(ocp.path.size())]) / pathScaling;
        pathMin = C.vertcat(
            [ocp.variable(pathVarNames[i]).min.getValue() for i in range(ocp.path.size())]) / pathScaling;
        pathConstraints.setOption("name", "PATH")
        pathConstraints.init()

        # pathConstraints.setInput(xOpt, 'x')
        pathConstraints.setInput([], 'z')
        # pathConstraints.setInput(uOpt, 'p')
        pathConstraints.setInput(0, 't')
        pathConstraints.evaluate()
        pathOpt = pathConstraints.getOutput()

        optimalValues = {}

        print 'min <= (name,optimal,nominal) <= max'
        for i in range(self.n_x):
            print ocp.variable(ocp.x[i].getName()).min.getValue(), ' <= (', ocp.x[i].getName(), ',', xOpt[i] * \
                                                                                                     stateScaling[
                                                                                                         i], ',', \
            stateScaling[i], ') <= ', ocp.variable(ocp.x[i].getName()).max.getValue()
            optimalValues[ocp.x[i].getName()] = xOpt[i] * stateScaling[i]

        for i in range(self.n_u):
            print ocp.variable(ocp.u[i].getName()).min.getValue(), ' <= (', ocp.u[i].getName(), ',', uOpt[i] * \
                                                                                                     controlScaling[
                                                                                                         i], ',', \
            controlScaling[i], ')  <= ', ocp.variable(ocp.u[i].getName()).max.getValue()
            if addDeltaUCons:
                print -deltaUCons[i] * controlScaling[i], ' <= (Delta(', ocp.u[i].getName(), ')/DTMPC,', 0, ',', \
                controlScaling[i], ')  <= ', deltaUCons[i] * controlScaling[i]

            optimalValues[ocp.u[i].getName()] = uOpt[i] * controlScaling[i]

        for i in range(len(pathVarNames)):
            print ocp.variable(pathVarNames[i]).min.getValue(), ' <= (', pathVarNames[i], ',', pathOpt[i] * pathScaling[
                i], ',', pathScaling[i], ') <= ', ocp.variable(pathVarNames[i]).max.getValue()
            optimalValues[pathVarNames[i]] = pathOpt[i] * pathScaling[i]

        plotTags = [ocp.x[i].getName() for i in range(ocp.x.size())]
        plotTags = plotTags + [ocp.u[i].getName() for i in range(ocp.u.size())]
        plotTags = plotTags + [sv.getName() for sv in ocp.beq(ocp.path)]

        self.plotTags = plotTags
        self.optimalValues = optimalValues

        # Constraint functions
        g = []
        g_min = []
        g_max = []

        self.XU0 = C.MX.sym("XU0", self.n_x + self.n_u, 1)
        Z = self.XU0[0:self.n_x]

        U0 = self.XU0[self.n_x:self.n_x + self.n_u]

        # Build up a graph of integrator calls
        obj = 0
        zf = C.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())]) / algStateScaling
        for k in range(self.n_k):
            Z, QF, zf = C.integratorOut(G(C.integratorIn(x0=Z, p=self.U[k], z0=zf)), "xf", "qf", "zf")

            errU = self.U[k] - U0
            obj = obj + QF #+ C.mul(C.mul(errU.T, controlCost), errU)

            U0 = self.U[k]

            [pathCons] = pathConstraints.call(C.daeIn(t=[], x=Z, z=zf, p=self.U[k]));
            g.append(pathCons)  ## be carefull on giving all inputs
            g_max.append(pathMax)
            g_min.append(pathMin)
            if addDeltaUCons:
                g.append(errU)
                g_max.append(deltaUCons * DT)
                g_min.append(-deltaUCons * DT)

        # errU = (self.U[-1]-uOpt)
        # errX = Z-xOpt
        # obj = obj + finalStateCost*C.mul((errX).trans(),(errX)) + C.mul(C.mul(errU.T,controlCost),errU)
        self.obj = obj

        ### Constrains
        g = C.vertcat(g)

        nlp = C.MXFunction(C.nlpIn(x=self.V, p=self.XU0), C.nlpOut(f=obj, g=g))
        nlp.init()
        self.odeF = odeF
        self.nlp = nlp

        solver = C.NlpSolver('ipopt', nlp)

        # remove the comment to implement the hessian
        solver.setOption('hessian_approximation', 'limited-memory')  # comment for exact hessian
        solver.setOption('print_user_options', 'yes')

        solver.setOption("tol", self.IPOPT_tol)  # IPOPT tolerance
        solver.setOption("dual_inf_tol", self.IPOPT_dual_inf_tol)  # dual infeasibility
        solver.setOption("constr_viol_tol", self.IPOPT_constr_viol_tol)  # primal infeasibility
        solver.setOption("compl_inf_tol", self.IPOPT_compl_inf_tol)  # complementarity
        #        solver.setOption("acceptable_tol",0.01)
        #        solver.setOption("acceptable_obj_change_tol",1e-6)
        #        solver.setOption("acceptable_constr_viol_tol",1e-6)
        solver.setOption("print_level", self.IPOPT_print_level)
        solver.setOption("max_iter", self.IPOPT_max_iter)  # IPOPT maximum iterations
        solver.setOption("max_cpu_time", self.IPOPT_max_cpu_time)  # IPOPT maximum iterations

        solver.init()

        ### Variable Bounds and initial guess
        solver.setInput(UMIN, 'lbx')  # u_L
        solver.setInput(UMAX, 'ubx')  # u_U
        solver.setInput(C.vertcat(g_min), 'lbg')  # g_L
        solver.setInput(C.vertcat(g_max), 'ubg')  # g_U

        self.solver = solver;

        u0N = C.vertcat([self.ocp.variable(self.ocp.u[i].getName()).initialGuess.getValue() for i in
                         range(self.n_u)]) / controlScaling

        USOL = [u0N for i in range(self.n_k)]

        self.USOL = USOL

    def solveProblem(self, x0, u0):

        # x0 = self.modifyStates(x0)

        x0 = x0 / self.stateScaling
        u0 = u0 / self.controlScaling

        V0 = self.mergeVariables(self.USOL)

        self.solver.setInput(V0, 'x0')
        self.solver.setInput(C.vertcat([x0, u0]), 'p')

        self.solver.evaluate()  # Solve NLP

        # =============================================================================
        #         if self.solver.getStat('iter_count') <= 1 or self.solver.getStat('return_status') == 'Invalid_Number_Detected':
        #             plotTags = self.plotTags
        #             try:
        #                 self.plotNLP(x0*self.stateScaling,plotTags,DTplot=1)
        #             except:
        #                 print "Unexpected error:", sys.exc_info()[0]
        #
        #             self.plotFigures()
        #
        #             Config = ConfigParser.ConfigParser()
        #             Config.read('config.ini')
        #             if Config.getboolean('MultipleShooting','plotonlyonce'):
        #                 Config.set('MultipleShooting','plotsolution',False)
        #                 cfgfile = open('config.ini','w')
        #                 Config.write(cfgfile)
        #                 cfgfile.close()
        #
        #             sameWindowPlots = Config.get('MultipleShooting','plotSameWindow')
        #             sameWindowPlots = sameWindowPlots.split(',')
        #
        #             for pw in sameWindowPlots:
        #                 tagsSamePlot = Config.get('MultipleShooting',pw)
        #                 tagsSamePlot = tagsSamePlot.split(',')
        #                 self.plotFiguresSameWindow(tagsSamePlot)
        #
        #             try:
        #                 msg = raw_input('Something went wrong, press enter to continue anyway...')
        #             except:
        #                 print "Unexpected error:", sys.exc_info()[0]
        # =============================================================================

        # Obtain results from optimization
        vSol = self.solver.output('x')
        objValue = self.solver.output('f')
        stats = self.solver.getStats()

        self.USOL = self.splitVariables(vSol)

        Config = ConfigParser.ConfigParser()
        Config.read('config.ini')
        if Config.getboolean('MultipleShooting', 'plotsolution'):
            plotTags = self.plotTags
            try:
                self.plotNLP(x0 * self.stateScaling, DTplot=1)
            except:
                print "Unexpected error:", sys.exc_info()[0]

            self.plotFigures()

            if Config.getboolean('MultipleShooting', 'plotonlyonce'):
                Config.set('MultipleShooting', 'plotsolution', False)
                cfgfile = open('config.ini', 'w')
                Config.write(cfgfile)
                cfgfile.close()

            sameWindowPlots = Config.get('MultipleShooting', 'plotSameWindow')
            sameWindowPlots = sameWindowPlots.split(',')

            for pw in sameWindowPlots:
                tagsSamePlot = Config.get('MultipleShooting', pw)
                tagsSamePlot = tagsSamePlot.split(',')
                self.plotFiguresSameWindow(tagsSamePlot, name=pw)

            if Config.getboolean('MultipleShooting', 'waitIfOptimized'):
                try:
                    msg = raw_input('The optimizer output has been plot, press enter to continue anyway...')
                except:
                    pass

        return self.USOL[0] * self.controlScaling, objValue, stats;

    def splitVariables(self, V):
        U = [V[k * self.n_u:(k + 1) * self.n_u] for k in range(self.n_k)]

        return U

    def mergeVariables(self, U):
        V = C.vertcat(U)
        return V

    def moveHorizion(self):

        self.USOL.append(NP.copy(self.USOL[-1]))
        del self.USOL[0]

    def plotNLP(self, x0, DTplot=10, additionalplottags=None):

        plotH = []
        plotAdd = []
        plotT = []
        plotTags = self.plotTags

        ocp = self.ocp
        algStateScaling = self.algStateScaling
        controlScaling = self.controlScaling
        stateScaling = self.stateScaling
        DT = self.DT
        Z = x0 / (self.stateScaling)
        odeF = self.odeF
        xOpt = self.xOpt
        uOpt = self.uOpt

        C.Integrator.loadPlugin("idas")
        G = C.Integrator("idas", odeF)
        G.setOption("reltol", self.INTG_REL_TOL)  # for CVODES and IDAS
        G.setOption("abstol", self.INTG_ABS_TOL)  # for CVODES and IDAS
        G.setOption("max_multistep_order", 5)  # for CVODES and IDAS
        G.setOption("max_step_size", self.IDAS_MAX_STEP_SIZE)  # for IDAS only
        G.setOption("tf", DTplot)
        G.init()

        pathIn = C.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t)

        if additionalplottags is None:
            Config = ConfigParser.ConfigParser()
            Config.read('config.ini')
            additionalplottags = Config.get('MultipleShooting', 'additionalplottags')
            additionalplottags = additionalplottags.split(',')

        self.additionalplottags = additionalplottags;

        addTagScale = C.vertcat([ocp.nominal(pv) for pv in additionalplottags])
        tagsFsx = C.substitute(C.vertcat([ocp.beq(sv) for sv in additionalplottags]), C.vertcat([ocp.x, ocp.z, ocp.u]),
                               C.vertcat([stateScaling * ocp.x, algStateScaling * ocp.z, controlScaling * ocp.u]))
        tagsF = C.SXFunction(pathIn, [tagsFsx]);
        tagsF.init()
        tagsF.setInput(xOpt, 'x')
        tagsF.setInput([], 'z')
        tagsF.setInput(uOpt, 'p')
        tagsF.setInput(0, 't')
        tagsF.evaluate()
        addTagsOpt = tagsF.getOutput()

        for i in range(len(additionalplottags)):
            self.optimalValues[additionalplottags[i]] = addTagsOpt[i]

        sxPlotTags = C.vertcat([ocp.beq(tag) for tag in plotTags])
        pathS = C.substitute(sxPlotTags, C.vertcat([ocp.x, ocp.z, ocp.u]),
                             C.vertcat([stateScaling * ocp.x, algStateScaling * ocp.z, controlScaling * ocp.u]))
        pathF = C.SXFunction(pathIn, [pathS, tagsFsx]);
        pathF.init()

        nPlot = NP.int(NP.round(self.DT / DTplot))
        zf = C.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())]) / algStateScaling
        for k in range(self.n_k):
            try:
                for jp in range(nPlot):
                    # Z,QF,zf = C.integratorOut(G(C.integratorIn(x0=Z,p=self.USOL[k],z0=zf)),"xf","qf","zf")
                    G.setInput(Z, 'x0')
                    G.setInput(self.USOL[k], 'p')
                    G.evaluate()
                    Z = G.getOutput('xf')
                    QF = G.getOutput('qf')
                    zf = G.getOutput('zf')
                    pathF.setInput(Z, 'x')
                    pathF.setInput(zf, 'z')
                    pathF.setInput(self.USOL[k], 'p')
                    t = k * DT + jp * DTplot
                    pathF.setInput(t, 't')
                    pathF.evaluate()
                    p = pathF.getOutput(0)
                    pAdd = pathF.getOutput(1)
                    plotH.append(p)
                    plotAdd.append(pAdd)
                    plotT.append(t)

            except:
                print 'something bad happened', sys.exc_info()[0]
                break

        plotDic = {}
        plotDic['t'] = plotT
        for si in range(len(plotTags)):
            plotDic[plotTags[si]] = [i[si] for i in plotH]
        for si in range(len(additionalplottags)):
            plotDic[additionalplottags[si]] = [i[si] for i in plotAdd]

        self.plotDic = plotDic

    #        self.plotFigures()

    def plotFigures(self):

        ocp = self.ocp
        plotDic = self.plotDic
        t = plotDic['t']
        plotTags = self.plotTags
        optimalValues = self.optimalValues;
        additionalplottags = self.additionalplottags

        plt.close('all')
        xMin = NP.zeros(len(plotTags))
        xMax = NP.zeros(len(plotTags))
        i = 0
        for si in plotTags:
            yiM = plotDic[si]

            minVal = ocp.variable(si).min.getValue()
            maxVal = ocp.variable(si).max.getValue()
            xMin[i] = minVal
            xMax[i] = maxVal
            i += 1

            optimalVal = optimalValues[si];

            tBounds = [t[0], t[-1]]
            minPlot = [minVal, minVal];
            maxPlot = [maxVal, maxVal];
            optimalPlot = [optimalVal, optimalVal]

            plt.figure()
            plt.plot(t, yiM, '-b', tBounds, minPlot, '--xr', tBounds, maxPlot, '-.or', tBounds, optimalPlot, ':Dk')

            # plt.legend([plotTags[si],'min','max','opt'])
            plt.title(si)
        print 'xMax:', xMax
        print 'xMin:', xMin
        for si in additionalplottags:
            yiM = plotDic[si]

            optimalVal = optimalValues[si];

            tBounds = [t[0], t[-1]]
            optimalPlot = [optimalVal, optimalVal]

            plt.figure()
            plt.plot(t, yiM, '-b', tBounds, optimalPlot, ':Dk')

            # plt.legend([plotTags[si],'min','max','opt'])
            plt.title(si)

        plt.ion()
        plt.show()

    def plotFiguresSameWindow(self, tags, name=''):

        plotDic = self.plotDic
        t = plotDic['t']

        plt.figure()

        for si in tags:
            plt.plot(t, plotDic[si])

        plt.legend(tags)
        plt.title(name)
        plt.ion()
        plt.show()

    def simulatePlot(self, x0, U=None, outputVars=None):

        if U is None:
            U = self.USOL

        outs = []
        tout = []

        if outputVars is not None:
            outputSX = C.vertcat([self.ocp.variable(outputVars[k]).beq for k in range(len(outputVars))])
            fout = C.SXFunction(C.daeIn(x=self.ocp.x, z=self.ocp.z, p=self.ocp.u, t=self.ocp.t), [outputSX])
            fout.init()

        zf = C.vertcat([self.ocp.variable(self.ocp.z[k].getName()).start for k in range(self.ocp.z.size())])

        xh = []
        zh = []

        xk = x0
        for k in range(len(U)):
            self.G.setInput(xk, 'x0')
            self.G.setInput(zf, 'z0')
            self.G.setInput(U[k], 'p')

            self.G.evaluate()

            xk = self.G.getOutput('xf')
            zf = self.G.getOutput('zf')

            xh.append(xk)
            zh.append(zf)

            tout.append((k + 1) * self.DT)
            if outputVars is not None:
                fout.setInput(xk, 'x')
                fout.setInput(zf, 'z')
                fout.setInput(U[k], 'p')
                fout.setInput((k + 1) * self.DT, 't')

                fout.evaluate()

                fouk = fout.getOutput()

                outs.append(fouk)

        return tout, xh, zh, outs

    def steadyStateOptimal(self, x0=None, u0=None, z0=None):

        ocp = self.ocp

        uMin = C.vertcat([ocp.variable(ocp.u[i].getName()).min.getValue() for i in range(ocp.u.size())])
        uMax = C.vertcat([ocp.variable(ocp.u[i].getName()).max.getValue() for i in range(ocp.u.size())])
        if u0 is None:
            u0 = C.vertcat([ocp.variable(ocp.u[i].getName()).initialGuess.getValue() for i in range(ocp.u.size())])

        xMin = C.vertcat([ocp.variable(ocp.x[i].getName()).min.getValue() for i in range(ocp.x.size())])
        xMax = C.vertcat([ocp.variable(ocp.x[i].getName()).max.getValue() for i in range(ocp.x.size())])
        if x0 is None:
            x0 = C.vertcat([ocp.variable(ocp.x[i].getName()).initialGuess.getValue() for i in range(ocp.x.size())])

        zMin = C.vertcat([ocp.variable(ocp.z[i].getName()).min.getValue() for i in range(ocp.z.size())])
        zMax = C.vertcat([ocp.variable(ocp.z[i].getName()).max.getValue() for i in range(ocp.z.size())])
        if z0 is None:
            z0 = C.vertcat([ocp.variable(ocp.z[i].getName()).initialGuess.getValue() for i in range(ocp.z.size())])

        ## PATH CONSTRAINTS
        pathIn = C.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t)
        pathConstraints = C.SXFunction(pathIn, [ocp.beq(ocp.path)]);
        pathMax = [ocp.variable(ocp.path[i].getName()).max.getValue() for i in range(ocp.path.size())];
        pathMin = [ocp.variable(ocp.path[i].getName()).min.getValue() for i in range(ocp.path.size())];
        pathConstraints.setOption("name", "PATH")
        pathConstraints.init()

        varList = [ocp.x, ocp.z, ocp.u]
        varMax = [xMax, zMax, uMax]
        varMin = [xMin, zMin, uMin]
        var0 = [x0, z0, u0]

        cons = [ocp.ode(ocp.x), ocp.alg, ocp.beq(ocp.path)]
        consMax = [NP.zeros(ocp.x.size()), NP.zeros(ocp.z.size()), C.vertcat(pathMax)]
        consMin = [NP.zeros(ocp.x.size()), NP.zeros(ocp.z.size()), C.vertcat(pathMin)]

        C.NlpSolver.loadPlugin('ipopt')
        nlp = C.SXFunction(C.nlpIn(x=C.vertcat(varList)), C.nlpOut(f=ocp.lterm, g=C.vertcat(cons)))
        nlp.init()

        solver = C.NlpSolver('ipopt', nlp)

        solver.setOption("max_iter", 100)  # IPOPT maximum iterations

        solver.init()

        ### Variable Bounds and initial guess
        solver.setInput(C.vertcat(varMin), 'lbx')
        solver.setInput(C.vertcat(varMax), 'ubx')

        solver.setInput(C.vertcat(consMin), 'lbg')  # g_L
        solver.setInput(C.vertcat(consMax), 'ubg')  # g_U

        solver.setInput(C.vertcat(var0), 'x0')  # g_U

        solver.evaluate()

        outX = solver.getOutput()

        splitF = C.SXFunction([C.vertcat(varList)], varList)
        splitF.init()
        splitF.setInput(outX)
        splitF.evaluate()

        xOpt = splitF.getOutput(0)
        zOpt = splitF.getOutput(1)
        uOpt = splitF.getOutput(2)

        return uOpt, xOpt, zOpt

    def steadyState(self, u, T=0):

        ocp = self.ocp

        xMin = C.vertcat([ocp.variable(ocp.x[i].getName()).min.getValue() for i in range(ocp.x.size())])
        xMax = C.vertcat([ocp.variable(ocp.x[i].getName()).max.getValue() for i in range(ocp.x.size())])
        x0 = C.vertcat([ocp.variable(ocp.x[i].getName()).initialGuess.getValue() for i in range(ocp.x.size())])

        zMin = C.vertcat([ocp.variable(ocp.z[i].getName()).min.getValue() for i in range(ocp.z.size())])
        zMax = C.vertcat([ocp.variable(ocp.z[i].getName()).max.getValue() for i in range(ocp.z.size())])
        z0 = C.vertcat([ocp.variable(ocp.z[i].getName()).initialGuess.getValue() for i in range(ocp.z.size())])

        if T > 0:
            tout, xh, zh, outs = self.simulatePlot(x0, [u for k in range(T)], outputVars=None)
            x0 = xh[-1]
            z0 = zh[-1]

        varList = [ocp.x, ocp.z]
        varMax = [xMax, zMax]
        varMin = [xMin, zMin]
        var0 = [x0, z0]

        cons = [ocp.ode(ocp.x), ocp.alg]
        consMax = [NP.zeros(ocp.x.size()), NP.zeros(ocp.z.size())]
        consMin = [NP.zeros(ocp.x.size()), NP.zeros(ocp.z.size())]

        C.NlpSolver.loadPlugin('ipopt')
        nlp = C.SXFunction(C.nlpIn(x=C.vertcat(varList), p=ocp.u), C.nlpOut(f=C.SX(0), g=C.vertcat(cons)))
        nlp.init()

        solver = C.NlpSolver('ipopt', nlp)

        solver.setOption("max_iter", 100)  # IPOPT maximum iterations

        solver.init()

        ### Variable Bounds and initial guess
        solver.setInput(C.vertcat(varMin), 'lbx')
        solver.setInput(C.vertcat(varMax), 'ubx')

        solver.setInput(C.vertcat(consMin), 'lbg')  # g_L
        solver.setInput(C.vertcat(consMax), 'ubg')  # g_U

        solver.setInput(C.vertcat(var0), 'x0')
        solver.setInput(u, 'p')

        solver.evaluate()

        outX = solver.getOutput()

        splitF = C.SXFunction([C.vertcat(varList)], varList)
        splitF.init()
        splitF.setInput(outX)
        splitF.evaluate()

        xss = splitF.getOutput(0)
        zss = splitF.getOutput(1)

        return xss, zss

    def testNlpCalling(self, x0, u0):

        x0 = x0 / self.stateScaling
        u0 = u0 / self.controlScaling

        self.nlp.setInput(C.vertcat([x0, u0]), 'p')
        V0 = self.mergeVariables(self.USOL)
        self.nlp.setInput(V0, 'x')
        self.nlp.evaluate()
