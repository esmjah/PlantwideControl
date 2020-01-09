clc
clear all
close all

format shortG

load D3D4_FeedbackRTO_new.mat
yOlgaMeas_fbrto = yOlgaMeas;
yOlgaMoni_fbrto = yOlgaMoni;
yKfMoni_fbrto = yKfMoni;
u_fbrto = u_OnSOC;

%load SimData_NMPC_GOR_2018.08.25_Olga_fine.mat
load D3D4_NMPC_600_8.mat

J_NMPC_new = J_NMPC;
Ju_NMPC_new = Ju_NMPC;
simTime_new = simTime;
u_NMPC_new = u_NMPC;
yOlgaMeas_new = yOlgaMeas;
yOlgaMoni_new = yOlgaMoni;
yKfMoni_new = yKfMoni;

load SimData_NMPC_GOR_2018.08.25_Olga_fine.mat
%load D3D4_NMPC_1200_12.mat

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

J_ideal_nominal = -28.4330;
J_ideal_d3 = -28.8272;
J_ideal_d4 = -29.2391;

u1_ideal_nominal = 1.2963;
u1_ideal_d3 = 0.8467;
u1_ideal_d4 = 0.8455;

u2_ideal_nominal = 1.3245;
u2_ideal_d3 = 1.3146;
u2_ideal_d4 = 0.8302;


plotStates = 0;
plotInputs = 1;
plotOutputs = 1;
plotMonitors = 1;


nU = length(u_NMPC(1,:));

NIT = length(x(:,1));
dt = 100;
x_lim = [0 50];
Hrs = round(length(simTime)/36)-5;
N = Hrs*36 + 1;
N1 = 5*36+1;
N2 = N+5*36;
time = simTime(N1:N2)/3600-5;
J_ideal = [J_ideal_nominal*ones(1,5*36) linspace(J_ideal_nominal,J_ideal_d3,10*36) J_ideal_d3*ones(1,10*36) linspace(J_ideal_d3,J_ideal_d4,10*36) J_ideal_d4*ones(1,N-35*36)];
u1_ideal = [u1_ideal_nominal*ones(1,5*36) linspace(u1_ideal_nominal,u1_ideal_d3,10*36) u1_ideal_d3*ones(1,10*36) linspace(u1_ideal_d3,u1_ideal_d4,10*36) u1_ideal_d4*ones(1,N-35*36)];
u2_ideal = [u2_ideal_nominal*ones(1,5*36) linspace(u2_ideal_nominal,u2_ideal_d3,10*36) u2_ideal_d3*ones(1,10*36) linspace(u2_ideal_d3,u2_ideal_d4,10*36) u2_ideal_d4*ones(1,N-35*36)];

J_fbrto = J_OnSOC(N1:N2);
Loss_fbrto = J_fbrto-J_ideal;

J_NMPC2 = J_NMPC(N1:N2);
Loss_NMPC = J_NMPC2-J_ideal;

J_NMPC2_new = J_NMPC_new(N1:N2);
Loss_NMPC_new = J_NMPC2_new-J_ideal;


disp(['cost_NMPC_Nominal = ', num2str(J_NMPC2(5*36-1))])
disp(['Loss_NMPC_Nominal = ', num2str(Loss_NMPC(5*36-1))])
disp(['cost_NMPC_Nominal_new = ', num2str(J_NMPC2_new(5*36-1))])
disp(['Loss_NMPC_Nominal_new = ', num2str(Loss_NMPC_new(5*36-1))])

disp(['cost_NMPC_d3 = ', num2str(J_NMPC2(25*36-1))])
disp(['Loss_NMPC_d3 = ', num2str(Loss_NMPC(25*36-1))])
disp(['cost_NMPC_d3_new = ', num2str(J_NMPC2_new(25*36-1))])
disp(['Loss_NMPC_d3_new = ', num2str(Loss_NMPC_new(25*36-1))])

disp(['cost_NMPC_d4 = ', num2str(J_NMPC2(50*36-1))])
disp(['Loss_NMPC_d4 = ', num2str(Loss_NMPC(50*36-1))])
disp(['cost_NMPC_d4_new = ', num2str(J_NMPC2_new(50*36-1))])
disp(['Loss_NMPC_d4_new = ', num2str(Loss_NMPC_new(50*36-1))])
%%
figure('Name','Cost')
plot(time,J_NMPC2,time,J_NMPC2_new,'-r',time,J_ideal,'--g')
legend('NMPC 1200','NMPC 600','Ideal')
%%
figure('Name','Loss')
clf
subplot(3,1,1)
plot(time,Loss_NMPC,time,Loss_NMPC_new,'r',time,zeros(size(time)),'--g','LineWidth',1)
title('Loss')
xlim(x_lim)
ylim([-0.1 0.2])
legend('NMPC 1200','NMPC 600','ideal')
subplot(3,1,2)
plot(time,Ju_NMPC(N1:N2,1),time,Ju_NMPC_new(N1:N2,1),'r',time,zeros(size(time)),'--g','LineWidth',1)
title('Ju_1')
xlim(x_lim)
legend('NMPC 1200','NMPC 600','setpoint')
subplot(3,1,3)
plot(time,Ju_NMPC(N1:N2,2),time,Ju_NMPC_new(N1:N2,1),'r',time,zeros(size(time)),'--g','LineWidth',1)
title('Ju_2')
xlim(x_lim)
legend('NMPC 1200','NMPC 600','setpoint')
%%
figure('Name','Cost--Loss')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(time,J_NMPC2,'Color',col1,'LineWidth',1)
title('Cost','Interpreter','latex')
hold on
plot(time,J_NMPC2_new,'-','Color',col2,'LineWidth',1)
plot(time,J_fbrto,'-','Color',col3,'LineWidth',1)
plot(time,J_ideal,':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('EMPC (Interval=1200 s)','EMPC (Interval=600 s)','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','NorthEast','Interpreter','latex','NumColumns',1)
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
plot(time,Loss_NMPC_new,'-','Color',col2,'LineWidth',1)
plot(time,Loss_fbrto,'-','Color',col3,'LineWidth',1)
plot(time,zeros(size(time)),':','Color',[0 .5 0],'LineWidth',1)
leg2 = legend('EMPC (Interval=1200 s)','EMPC (Interval=600 s)','Feedback-RTO','Steady-State Optimal');
set(leg2,'Location','NorthWest','Interpreter','latex','NumColumns',1)
xlim(x_lim)
%ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('L [\$/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2_revised\Figures\cost_Loss_GOR
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\cost_Loss_GOR_comp
print -dpdf C:\Git\PlantwideControl\SavedResults\Figures\cost_Loss_GOR_comp
%%
figure('Name','Distirbances-GOR')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(time,yOlgaMoni_new(N1:N2,11),'-','Color',col4,'LineWidth',1)
title('Well A gas-oil ratio','Interpreter','latex')
hold on
plot(time,yKfMoni(N1:N2,11),'-','Color',col1,'LineWidth',1)
plot(time,yKfMoni_new(N1:N2,11),'-','Color',col2,'LineWidth',1)
plot(time,yKfMoni_fbrto(N1:N2,11),'-','Color',col3,'LineWidth',1)
leg1 = legend('Olga Simulator','EMPC (Interval=1200 s)','EMPC (Interval=600 s)','Feedback-RTO');
set(leg1,'Location','Best','Interpreter','latex')
xlim(x_lim)
ylim([-0.005 0.035])
xlabel('time [h]','Interpreter','latex')
ylabel('$GOR_{\textrm{A}}$ [-]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
subplot(2,1,2,'OuterPosition',[0 0 1 1/2])
plot(time,yOlgaMoni_new(N1:N2,12),'-','Color',col4,'LineWidth',1)
title('Well B gas-oil ratio','Interpreter','latex')
hold on
plot(time,yKfMoni(N1:N2,12),'-','Color',col1,'LineWidth',1)
plot(time,yKfMoni_new(N1:N2,12),'-','Color',col2,'LineWidth',1)
plot(time,yKfMoni_fbrto(N1:N2,12),'-','Color',col3,'LineWidth',1)
leg2 = legend('Olga Simulator','EMPC (Interval=1200 s)','EMPC (Interval=600 s)','Feedback-RTO');
set(leg2,'Location','Best','Interpreter','latex')
xlim(x_lim)
ylim([-0.005 0.035])
xlabel('time [h]','Interpreter','latex')
ylabel('$GOR_{\textrm{B}}$ [-]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';
%print -depsc Z:\Dropbox\OptimumSeeking\Manuscript2_revised\Figures\disturbances_GOR
print -depsc C:\Git\PlantwideControl\SavedResults\Figures\disturbances_GOR_comp
print -dpdf C:\Git\PlantwideControl\SavedResults\Figures\disturbances_GOR_coimp
%%
figure('Name','Inputs--Gor')
clf
rect = [0, 0, 12, 10];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 10],'PaperPosition',rect)
subplot(2,1,1,'OuterPosition',[0 1/2 1 1/2])
plot(time,u_NMPC(N1:N2,1),'Color',col1,'LineWidth',1)
title('Well A gas injection rate','Interpreter','latex')
hold on
plot(time,u_NMPC_new(N1:N2,1),'-','Color',col2,'LineWidth',1)
plot(time,u_fbrto(N1:N2,1),'-','Color',col3,'LineWidth',1)
plot(time,u1_ideal,':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('EMPC (Interval=1200 s)','EMPC (Interval=600 s)','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','NorthEast','Interpreter','latex','NumColumns',1)
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
plot(time,u_NMPC_new(N1:N2,2),'-','Color',col2,'LineWidth',1)
plot(time,u_fbrto(N1:N2,2),'-','Color',col3,'LineWidth',1)
plot(time,u2_ideal,':','Color',[0 .5 0],'LineWidth',1)
leg1 = legend('EMPC (Interval=1200 s)','EMPC (Interval=600 s)','Feedback-RTO','Steady-State Optimal');
set(leg1,'Location','SouthWest','Interpreter','latex','NumColumns',1)
xlim(x_lim)
%ylim([0.9 1.4])
xlabel('time [h]','Interpreter','latex')
ylabel('$w_{\textrm{inj,B}}$ [kg/s]','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';

print -depsc C:\Git\PlantwideControl\SavedResults\Figures\injections_GOR_comp
print -dpdf C:\Git\PlantwideControl\SavedResults\Figures\injections_GOR_comp
%%
figure('Name','NMPC Pressure Controls')
clf
rect = [0, 0, 12, 14];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 14],'PaperPosition',rect)

subplot(3,1,1,'OuterPosition',[0 2/3 1 1/3])
plot(time,u_NMPC(N1:N2,3),'--','Color',[0 .5 0],'LineWidth',1)
hold on
plot(time,yOlgaMeas(N1:N2,9),'-','Color',col1,'LineWidth',1)
plot(time,yOlgaMeas_new(N1:N2,9),'-','Color',col2,'LineWidth',1)
plot(time,yOlgaMeas_fbrto(N1:N2,9),'-','Color',col3,'LineWidth',1)
title('Pressure drop of topside valve','Interpreter','latex')
xlim(x_lim)
ylim([5e4 9e4])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{top}}$ [Pa]','Interpreter','latex')
leg1 = legend('Setpoint','Interval=1200 s','Interval=600 s', 'Feedback-RTO');
set(leg1,'Location','SouthEast','Interpreter','latex','NumColumns',3)
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(3,1,2,'OuterPosition',[0 1/3 1 1/3])
plot(time,u_NMPC(N1:N2,4),'--','Color',[0 .5 0],'LineWidth',1)
hold on
plot(time,yOlgaMeas(N1:N2,7),'-','Color',col1,'LineWidth',1)
plot(time,yOlgaMeas_new(N1:N2,7),'-','Color',col2,'LineWidth',1)
plot(time,yOlgaMeas_fbrto(N1:N2,7),'-','Color',col3,'LineWidth',1)
title('Pressure drop of well A wellhead valve','Interpreter','latex')
xlim(x_lim)
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,A}}$ [Pa]','Interpreter','latex')
leg2 = legend('Setpoint','Interval=1200 s','Interval=600 s', 'Feedback-RTO');
set(leg2,'Location','Best','Interpreter','latex','NumColumns',3)
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(3,1,3,'OuterPosition',[0 0 1 1/3])
plot(time,u_NMPC(N1:N2,5),'--','Color',[0 .5 0],'LineWidth',1)
hold on
plot(time,yOlgaMeas(N1:N2,8),'-','Color',col1,'LineWidth',1)
plot(time,yOlgaMeas_new(N1:N2,8),'-','Color',col2,'LineWidth',1)
plot(time,yOlgaMeas_fbrto(N1:N2,8),'-','Color',col3,'LineWidth',1)
title('Pressure drop of well B wellhead valve','Interpreter','latex')
xlim(x_lim)
ylim([1.8e5 2.5e5])
xlabel('time [h]','Interpreter','latex')
ylabel('$\Delta P_{\textrm{wh,B}}$ [Pa]','Interpreter','latex')
leg3 = legend('Setpoint','Interval=1200 s','Interval=600 s', 'Feedback-RTO');
set(leg3,'Location','Best','Interpreter','latex','NumColumns',1)
axs = gca;
axs.TickLabelInterpreter = 'latex';

print -depsc C:\Git\PlantwideControl\SavedResults\Figures\MPC_pressure_controls_comp
print -dpdf C:\Git\PlantwideControl\SavedResults\Figures\MPC_pressure_controls_cmp
