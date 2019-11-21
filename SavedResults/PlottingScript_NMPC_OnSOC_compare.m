clc
clear all
close all
%load('SimData_NCO_d1_2017.05.18_testModel.mat')
load SimData_SOC2_Nominal_2019.03.02.mat
simTime_SOC = simTime/3600;

load SimData_OnSOC_Nominal_2018.06.25_Olga.mat
load SimData_NMPC_Nominal_2018.06.24_Olga.mat


plotStates = 0;
plotInputs = 1;
plotOutputs = 0;
plotMonitors = 0;

J_opt = -28.433;

col1 = [8, 118, 184]/255;
col2 = [36, 162, 14]/255;
col3 = [233, 11, 0]/255;

figure('Name','getDefaultColors')
h = plot(eye(6));
c = get(h,'Color');
close getDefaultColors

col1 = c{1};
col2 = c{2};
col3 = c{3};
col4 = c{4};
col5 = c{5};
col6 = c{6};
col6 = [0 0 0];

u = u_NMPC;

nU = length(u(1,:));

NIT = length(x(:,1));
dt = 100;
time = simTime/3600;%(0:dt:(NIT-1)*dt)'/3600;
x_lim = [0 30];
N = length(time);

if(plotStates)
    nX = length(x(1,:));
    for ix=1:nX
        figure('Name',states{ix})
        clf
        plot(time,x(:,ix),'LineWidth',1)
        xlim(x_lim)
        legend(states{ix})
        hold on
        plot(time,xMin(ix)*ones(NIT,1),'-r','LineWidth',2)
        xlim(x_lim)
        plot(time,xMax(ix)*ones(NIT,1),'-r','LineWidth',2)
        plot(time,xOptimal(ix)*ones(NIT,1),'--','Color',[0 .5 0],'LineWidth',2)
        xlim(x_lim)
    end
end

%%
figure('Name','Cost')
clf
rect = [0, 0, 12, 14];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 14],'PaperPosition',rect)
subplot(3,1,1,'OuterPosition',[0 2/3 1 1/3])
plot(time(1:N),J_NMPC(1:N),'Color',col1,'LineWidth',1)
hold on
plot(simTime_SOC,J_SOC,'-','Color',col3,'LineWidth',1)
plot(time(1:N),J_OnSOC(1:N),'--','Color',col2,'LineWidth',1)
plot(time,J_opt*ones(NIT,1),':','Color',[0 .5 0],'LineWidth',1)
title('Cost','Interpreter','latex')
xlim(x_lim)
xlabel('time [h]','Interpreter','latex')
ylabel('J [\$/s]','Interpreter','latex')
leg1 = legend('EMPC','Self-Optimizing','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','NorthEast','Interpreter','latex','NumColumns',2)
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(3,1,2,'OuterPosition',[0 1/3 1 1/3])
plot(time(1:N),Ju_NMPC(1:N,1),'-','Color',col1,'LineWidth',1)
hold on
plot(simTime_SOC,Ju_SOC(:,1),'-','Color',col3,'LineWidth',1)
plot(time(1:N),Ju_OnSOC(1:N,1),'--','Color',col2,'LineWidth',1)
plot(time,zeros(NIT,1),':','Color',[0 .5 0],'LineWidth',1)
title('Gradient with respect to well A injection rate','Interpreter','latex')
xlim(x_lim)
xlabel('time [h]','Interpreter','latex')
ylabel('$J_{\textrm{u,A}} $[\$/kg]','Interpreter','latex')
leg1 = legend('EMPC','Self-Optimizing','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','SouthEast','Interpreter','latex','NumColumns',2)
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(3,1,3,'OuterPosition',[0 0 1 1/3])
plot(time(1:N),Ju_NMPC(1:N,2),'-','Color',col1,'LineWidth',1)
hold on
plot(simTime_SOC,Ju_SOC(:,2),'-','Color',col3,'LineWidth',1)
plot(time(1:N),Ju_OnSOC(1:N,2),'--','Color',col2,'LineWidth',1)
plot(time,zeros(NIT,1),':','Color',[0 .5 0],'LineWidth',1)
title('Gradient with respect to well B injection rate','Interpreter','latex')
xlim(x_lim)
xlabel('time [h]','Interpreter','latex')
ylabel('$J_{\textrm{u,B}} $[\$/kg]','Interpreter','latex')
leg1 = legend('EMPC','Self-Optimizing','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','SouthEast','Interpreter','latex','NumColumns',2)
axs = gca;
axs.TickLabelInterpreter = 'latex';

%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2\Figures\nominal_cost
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\nominal_cost
%%
figure('Name','DOF')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(time(1:N),u_NMPC(1:N,1),'-','Color',col1,'LineWidth',1)
title('Well A gas injection rate','Interpreter','latex')
hold on
plot(simTime_SOC,u_SOC(:,1),'Color',col3,'LineWidth',1)
plot(time(1:N),u_OnSOC(1:N,1),'--','Color',col2,'LineWidth',1)
plot(time,uOptimal(1)*ones(NIT,1),':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('EMPC','Self-Optimizing','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','SouthEast','Interpreter','latex','NumColumns',2)
xlim(x_lim)
ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('$w_{\textrm{inj,A}}$ [kg/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
subplot(2,1,2,'OuterPosition',[0 0 1 1/2])
plot(time(1:N),u_NMPC(1:N,2),'-','Color',col1,'LineWidth',1)
title('Well B gas injection rate','Interpreter','latex')
hold on
plot(simTime_SOC,u_SOC(:,2),'Color',col3,'LineWidth',1)
plot(time(1:N),u_OnSOC(1:N,2),'--','Color',col2,'LineWidth',1)
plot(time,uOptimal(2)*ones(NIT,1),':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('EMPC','Self-Optimizing','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','SouthEast','Interpreter','latex','NumColumns',2)
xlim(x_lim)
ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('$w_{\textrm{inj,B}}$ [kg/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2\Figures\nominal_injections
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\nominal_injections
%%
figure('Name','NMPC Pressure Controls')
clf
rect = [0, 0, 12, 14];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 14],'PaperPosition',rect)

subplot(3,1,1,'OuterPosition',[0 2/3 1 1/3])
plot(time(1:N),u_NMPC(1:N,3),'-',time(1:N),yOlgaMeas(:,9),'-','LineWidth',1)
hold on
plot(time,yKfMeas(:,9),'--','Color',[0 .5 0],'LineWidth',1)
plot(time,uOptimal(3)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of topside valve','Interpreter','latex')
xlim(x_lim)
ylim([5e4 9e4])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{top}}$ [Pa]','Interpreter','latex')
leg1 = legend('Optimal Setpoint','Olga Simulator','Dynamic Model','Active Constraint');
set(leg1,'Location','NorthEast','Interpreter','latex','NumColumns',2)
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(3,1,2,'OuterPosition',[0 1/3 1 1/3])
plot(time(1:N),u_NMPC(1:N,4),'-',time(1:N),yOlgaMeas(:,7),'-','LineWidth',1)
hold on
plot(time,yKfMeas(:,7),'--','Color',[0 .5 0],'LineWidth',1)
plot(time,uOptimal(4)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of well A wellhead valve','Interpreter','latex')
xlim(x_lim)
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,A}}$ [Pa]','Interpreter','latex')
leg2 = legend('Optimal Setpoint','Olga Simulator','Dynamic Model','Active Constraint');
set(leg2,'Location','NorthEast','Interpreter','latex','NumColumns',2)
subplot(3,1,3,'OuterPosition',[0 0 1 1/3])
axs = gca;
axs.TickLabelInterpreter = 'latex';

plot(time(1:N),u_NMPC(1:N,5),'-',time(1:N),yOlgaMeas(:,8),'-','LineWidth',1)
hold on
plot(time,yKfMeas(:,8),'--','Color',[0 .5 0],'LineWidth',1)
plot(time,uOptimal(5)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of well B wellhead valve','Interpreter','latex')
xlim(x_lim)
ylim([1.8e5 3e5])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,B}}$ [Pa]','Interpreter','latex')
leg3 = legend('Optimal Setpoint','Olga Simulator','Dynamic Model','Active Constraint');
set(leg3,'Location','NorthEast','Interpreter','latex','NumColumns',2)
axs = gca;
axs.TickLabelInterpreter = 'latex';

%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2\Figures\constrained_inputs
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\MPC_pressure_controls
%%
load SimData_OnSOC_Nominal_2018.06.25_Olga.mat
figure('Name','Feedback RTO Pressure Controls')
clf
rect = [0, 0, 12, 14];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 14],'PaperPosition',rect)

subplot(3,1,1,'OuterPosition',[0 2/3 1 1/3])
plot(time(1:N),u_OnSOC(1:N,3),'-',time(1:N),yOlgaMeas(1:N,9),'-','LineWidth',1)
hold on
plot(time,yKfMeas(1:N,9),'--','Color',[0 .5 0],'LineWidth',1)
plot(time,uOptimal(3)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of topside valve','Interpreter','latex')
xlim(x_lim)
ylim([5e4 9e4])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{top}}$ [Pa]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
leg1 = legend('Optimal Setpoint','Olga Simulator','Dynamic Model','Active Constraint');
set(leg1,'Location','NorthEast','Interpreter','latex','NumColumns',2)

subplot(3,1,2,'OuterPosition',[0 1/3 1 1/3])
plot(time(1:N),u_OnSOC(1:N,4),'-',time(1:N),yOlgaMeas(1:N,7),'-','LineWidth',1)
hold on
plot(time,yKfMeas(1:N,7),'--','Color',[0 .5 0],'LineWidth',1)
plot(time,uOptimal(4)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of well A wellhead valve','Interpreter','latex')
xlim(x_lim)
ylim([1.8e5 3e5])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,A}}$ [Pa]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
leg2 = legend('Optimal Setpoint','Olga Simulator','Dynamic Model','Active Constraint');
set(leg2,'Location','NorthEast','Interpreter','latex','NumColumns',2)

subplot(3,1,3,'OuterPosition',[0 0 1 1/3])
plot(time(1:N),u_OnSOC(1:N,5),'-',time(1:N),yOlgaMeas(1:N,8),'-','LineWidth',1)
hold on
plot(time,yKfMeas(1:N,8),'--','Color',[0 .5 0],'LineWidth',1)
plot(time,uOptimal(5)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of well B wellhead valve','Interpreter','latex')
xlim(x_lim)
ylim([1.8e5 3.4e5])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,B}}$ [Pa]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
leg3 = legend('Optimal Setpoint','Olga Simulator','Dynamic Model','Active Constraint');
set(leg3,'Location','NorthEast','Interpreter','latex','NumColumns',2)
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2\Figures\constrained_inputs
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\FeedbackRTO_pressure_controls

%%
load SimData_SOC2_Nominal_2019.03.02.mat
simTime_SOC = simTime/3600;
figure('Name','Self-Optimizing Pressure Controls')
clf
rect = [0, 0, 12, 14];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 14],'PaperPosition',rect)

subplot(3,1,1,'OuterPosition',[0 2/3 1 1/3])
plot(simTime_SOC,u_SOC(:,3),'-',simTime_SOC,yOlgaMeas(:,9),'-','LineWidth',1)
hold on
plot(time,uOptimal(3)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of topside valve','Interpreter','latex')
xlim(x_lim)
ylim([5e4 9e4])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{top}}$ [Pa]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
leg1 = legend('Optimal Setpoint','Olga Simulator','Active Constraint');
set(leg1,'Location','NorthEast','Interpreter','latex')

subplot(3,1,2,'OuterPosition',[0 1/3 1 1/3])
plot(simTime_SOC,u_SOC(:,4),'-',simTime_SOC,yOlgaMeas(:,7),'-','LineWidth',1)
hold on
plot(time,uOptimal(4)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of well A wellhead valve','Interpreter','latex')
xlim(x_lim)
ylim([1.8e5 3e5])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,A}}$ [Pa]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
leg2 = legend('Optimal Setpoint','Olga Simulator','Active Constraint');
set(leg2,'Location','NorthEast','Interpreter','latex')

subplot(3,1,3,'OuterPosition',[0 0 1 1/3])
plot(simTime_SOC,u_SOC(:,5),'-',simTime_SOC,yOlgaMeas(:,8),'-','LineWidth',1)
hold on
plot(time,uOptimal(5)*ones(NIT,1),':k','LineWidth',1)
title('Pressure drop of well B wellhead valve','Interpreter','latex')
xlim(x_lim)
ylim([1.8e5 3.4e5])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,B}}$ [Pa]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
leg3 = legend('Optimal Setpoint','Olga Simulator','Active Constraint');
set(leg3,'Location','NorthEast','Interpreter','latex')
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\SOC_pressure_controls
%%
figure('Name','SOC_CVs')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(simTime_SOC,y1,'-','Color',col3,'LineWidth',1)
title('Gas flow rate at wellhead well A','Interpreter','latex')
hold on
plot(time,r1*ones(NIT,1),':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('Controlled Variable','Optimal Setpoint');
set(leg1,'Location','SouthEast','Interpreter','latex')
xlim(x_lim)
ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('$w_{\textrm{Gwh,A}}$ [kg/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
subplot(2,1,2,'OuterPosition',[0 0 1 1/2])
plot(simTime_SOC,y2,'-','Color',col3,'LineWidth',1)
title('Gas flow rate at wellhead well B','Interpreter','latex')
hold on
plot(time,r2*ones(NIT,1),':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('Controlled Variable','Optimal Setpoint');
set(leg1,'Location','SouthEast','Interpreter','latex')
xlim(x_lim)
ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('$w_{\textrm{Gwh,B}}$ [kg/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2\Figures\nominal_injections
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\nominal_SOC_CVs
%%
plotInputs = 0;
if(plotInputs)
    nu = length(u(1,:));
    for iu=1:nu
        figure('Name',controlTags{iu})
        clf
        plot(time(1:N),u_NMPC(1:N,iu),time(1:N),u_OnSOC(1:N,iu),'--','LineWidth',1)
        title(controlTags{iu})
        legend('EMPC','Dynamic Feedback-RTO')
        hold on
        %plot(time,uMin(iu)*ones(NIT,1),'-r','LineWidth',2)
        %plot(time,uMax(iu)*ones(NIT,1),'-r','LineWidth',2)
        plot(time,uOptimal(iu)*ones(NIT,1),'--','Color',[0 .5 0],'LineWidth',2)
        xlim(x_lim)
    end
end
plotOutputs = 0;
if(plotOutputs)
    ny = length(yOlgaMeas(1,:));
    idx = [2 3 1]+nu-3;
    for iy=1:ny
        figure('Name',measurementTags{iy})
        clf
        plot(time,yOlgaMeas(:,iy),'LineWidth',1)
        hold on
        plot(time,yKfMeas(:,iy),'--','Color',col3,'LineWidth',1)
        legend(measurementTags{iy},measurementTagsKf{iy})
        if(iy>(ny-3))
            plot(time,u(:,idx(iy-(ny-3))),'--k','LineWidth',1)
        end
        xlim(x_lim)
    end
end
plotMonitors = 0;
if(plotMonitors)
    ny = length(yOlgaMoni(1,:));
    for iy=1:ny
        figure('Name',measurementTagsMonitor{iy})
        clf
        plot(time,yOlgaMoni(:,iy),'LineWidth',1)
        hold on
        plot(time,yKfMoni(:,iy),'--','Color',col2,'LineWidth',1)
        legend(measurementTagsMonitor{iy},measurementTagsMonitorKf{iy})
    end
    xlim(x_lim)
end