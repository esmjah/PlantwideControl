clc
clear all
close all

format shortG

% Pres distirbance
load D1D2_SelfOptimizing.mat
load D1D2_FeedbackRTO_new.mat
load D1D2_NMPC_1200_12.mat

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

J_ideal_nominal = -28.433;
J_ideal_d1 = -27.8784;
J_ideal_d2 = -27.3130;

u1_ideal_nominal = 1.2963;
u1_ideal_d3 = 1.2788;
u1_ideal_d4 = 1.2765;

u2_ideal_nominal = 1.3245;
u2_ideal_d3 = 1.3209;
u2_ideal_d4 = 1.3071;

plotStates = 0;
plotInputs = 1;
plotOutputs = 1;
plotMonitors = 1;


nU = length(u_SOC(1,:));

NIT = length(x(:,1));
dt = 100;
x_lim = [0 50];
Hrs = round(length(simTime)/36)-5;
N = Hrs*36 + 1;
N1 = 5*36+1;
N2 = N+5*36;
time = simTime(N1:N2)/3600-5;
J_ideal = [J_ideal_nominal*ones(1,5*36) linspace(J_ideal_nominal,J_ideal_d1,10*36) J_ideal_d1*ones(1,10*36) linspace(J_ideal_d1,J_ideal_d2,10*36) J_ideal_d2*ones(1,N-35*36)];
u1_ideal = [u1_ideal_nominal*ones(1,5*36) linspace(u1_ideal_nominal,u1_ideal_d3,10*36) u1_ideal_d3*ones(1,10*36) linspace(u1_ideal_d3,u1_ideal_d4,10*36) u1_ideal_d4*ones(1,N-35*36)];
u2_ideal = [u2_ideal_nominal*ones(1,5*36) linspace(u2_ideal_nominal,u2_ideal_d3,10*36) u2_ideal_d3*ones(1,10*36) linspace(u2_ideal_d3,u2_ideal_d4,10*36) u2_ideal_d4*ones(1,N-35*36)];

J_SOC2 = J_SOC(N1:N2);
J_OnSOC2 = J_OnSOC(N1:N2);
J_NMPC2 = J_NMPC(N1:N2);
Loss_NMPC = J_NMPC2-J_ideal;
Loss_OnSOC = J_OnSOC2-J_ideal;
Loss_SOC = J_SOC2-J_ideal;

Loss_NMPC_Nominal = Loss_NMPC(5*36-1)
Loss_OnSOC_Nominal = Loss_OnSOC(5*36-1)
Loss_SOC_Nominal = Loss_SOC(5*36-1)

Loss_NMPC_d1 = Loss_NMPC(25*36-1)
Loss_OnSOC_d1 = Loss_OnSOC(25*36-1)
Loss_SOC_d1 = Loss_SOC(25*36-1)

Loss_NMPC_d2 = Loss_NMPC(50*36-1)
Loss_OnSOC_d2 = Loss_OnSOC(50*36-1)
Loss_SOC_d2 = Loss_SOC(50*36-1)
%%
figure('Name','Cost')
plot(time,J_OnSOC2,time,J_NMPC2,'--',time,J_SOC2,'--k',time,J_ideal,'--g')
legend('Proposed','NMPC','Self-Optimizng','Ideal')
%%
figure('Name','Loss')
clf
subplot(3,1,1)
plot(time,Loss_OnSOC,time,Loss_NMPC,'--',time,Loss_SOC,'--k',time,zeros(size(time)),'--g','LineWidth',1)
title('Loss')
xlim(x_lim)
ylim([-0.1 0.2])
legend('Proposed','NMPC','SOC','ideal')
subplot(3,1,2)
plot(time,Ju_OnSOC(N1:N2,1),time,Ju_NMPC(N1:N2,1),'--',time,Ju_SOC(N1:N2,1),'--k',time,zeros(size(time)),'--g','LineWidth',1)
title('Ju_1')
xlim(x_lim)
legend('Proposed','NMPC','SOC','setpoint')
subplot(3,1,3)
plot(time,Ju_OnSOC(N1:N2,2),time,Ju_NMPC(N1:N2,2),'--',time,Ju_SOC(N1:N2,1),'--k',time,zeros(size(time)),'--g','LineWidth',1)
title('Ju_2')
xlim(x_lim)
legend('Proposed','NMPC','SOC','setpoint')

%%
figure('Name','Cost--Loss')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(time,J_NMPC2,'Color',col1,'LineWidth',1)
title('Cost','Interpreter','latex')
hold on
plot(time,J_SOC2,'-','Color',col3,'LineWidth',1)
plot(time,J_OnSOC2,'--','Color',col2,'LineWidth',1)
plot(time,J_ideal,':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('EMPC','Self-Optimizng','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','Best','Interpreter','latex','NumColumns',2)
xlim(x_lim)
%ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('J [\$/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
subplot(2,1,2,'OuterPosition',[0 0 1 1/2])
plot(time,Loss_NMPC,'Color',col1,'LineWidth',1)
title('Loss','Interpreter','latex')
hold on
plot(time,Loss_SOC,'-','Color',col3,'LineWidth',1)
plot(time,Loss_OnSOC,'--','Color',col2,'LineWidth',1)
plot(time,zeros(size(time)),':','Color',[0 .5 0],'LineWidth',1)
leg2 = legend('EMPC','Self-Optimizng','Feedback-RTO','Steady-State Optimal');
set(leg2,'Location','Best','Interpreter','latex','NumColumns',2)
xlim(x_lim)
%ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('L [\$/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2_revised\Figures\cost_Loss_Pres
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\cost_Loss_Pres
print -dpdf C:\Git\PlantwideControl\SavedResults\Figures\cost_Loss_Pres
%%
figure('Name','Distirbances-P_res')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(time,yOlgaMoni(N1:N2,12),'LineWidth',1)
title('Well A reservoir pressure','Interpreter','latex')
hold on
plot(time,yKfMoni(N1:N2,12),'--','Color',col2,'LineWidth',1)
leg1 = legend('Actual (Olga Simulator)','EKF Estimation (Dynamic Model)');
set(leg1,'Location','Best','Interpreter','latex')
xlim(x_lim)
ylim([154 161])
xlabel('time [h]','Interpreter','latex')
ylabel('$P_{\textrm{res,A}}$ [bar]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
subplot(2,1,2,'OuterPosition',[0 0 1 1/2])
plot(time,yOlgaMoni(N1:N2,13),'LineWidth',1)
title('Well B reservoir pressure','Interpreter','latex')
hold on
plot(time,yKfMoni(N1:N2,13),'--','Color',col2,'LineWidth',1)
leg2 = legend('Actual (Olga Simulator)','EKF Estimation (Dynamic Model)');
set(leg2,'Location','Best','Interpreter','latex')
xlim(x_lim)
ylim([164 171])
xlabel('time [h]','Interpreter','latex')
ylabel('$P_{\textrm{res,B}}$ [bar]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2_revised\Figures\disturbances_Pres
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\disturbances_Pres
print -dpdf C:\Git\PlantwideControl\SavedResults\Figures\disturbances_Pres
%%
figure('Name','Inpus--Pres')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(time,u_NMPC(N1:N2,1),'Color',col1,'LineWidth',1)
title('Well A gas injection rate','Interpreter','latex')
hold on
plot(time,u_SOC(N1:N2,1),'-','Color',col3,'LineWidth',1)
plot(time,u_OnSOC(N1:N2,1),'--','Color',col2,'LineWidth',1)
plot(time,u1_ideal,':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend({'EMPC','Self-Optimizng','Feedback-RTO','Steady-State Optimal'});
set(leg1,'Location','Best','Interpreter','latex','NumColumns',2)
xlim(x_lim)
%ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('$w_{\textrm{inj,A}}$ [kg/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
subplot(2,1,2,'OuterPosition',[0 0 1 1/2])
plot(time,u_NMPC(N1:N2,2),'Color',col1,'LineWidth',1)
title('Well B gas injection rate','Interpreter','latex')
hold on
plot(time,u_SOC(N1:N2,2),'-','Color',col3,'LineWidth',1)
plot(time,u_OnSOC(N1:N2,2),'--','Color',col2,'LineWidth',1)
plot(time,u2_ideal,':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend({'EMPC','Self-Optimizng','Feedback-RTO','Steady-State Optimal'});
set(leg1,'Location','Best','Interpreter','latex','NumColumns',2)
xlim(x_lim)
%ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('$w_{\textrm{inj,B}}$ [kg/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2_revised\Figures\injections_Pres
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\injections_Pres
print -dpdf C:\Git\PlantwideControl\SavedResults\Figures\injections_Pres
%%
if(plotInputs)
    nu = length(u_NMPC(1,:));
    for iu=1:nu
        figure('Name',controlTags{iu})
        clf
        title(controlTags{iu})
        plot(time,u_NMPC(N1:N2,iu),time,u_OnSOC(N1:N2,iu),'--',time,u_SOC(N1:N2,iu),'--k','LineWidth',1)
        hold on
        %plot(time,uMin(iu)*ones(NIT,1),'-r','LineWidth',2)
        %plot(time,uMax(iu)*ones(NIT,1),'-r','LineWidth',2)
        plot(time,uOptimal(iu)*ones(size(time)),'--g','Color',[0 .5 0],'LineWidth',2)
        legend('NMPC','Proposed','SOC','Nominal')
        xlim(x_lim)
    end
end

figure('Name','SOC')
clf
subplot(2,1,1)
plot(time,y1(N1:N2),time,r1*ones(size(time)),'LineWidth',1)
legend('c1','r1')
subplot(2,1,2)
plot(time,y2(N1:N2),time,r2*ones(size(time)),'LineWidth',1)
legend('c2','r2')

