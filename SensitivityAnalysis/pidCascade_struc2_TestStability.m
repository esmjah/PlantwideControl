clc
clear all
close all
format long g
load jacobians_struc2_zOpt.mat
load('TuningValues.mat')
warning off

inputNames = {'inj1','inj2','zA1','zA2','zINL'}';
outputNames = {'yA1','yT1','yA2','yT2','yINL','yTOP'}';

sys = ss(A,B,C,D,'InputName',inputNames,'OutputName',outputNames);

s = tf('s');
%sys = sys_in;% = minreal(sys_in);

%pole(sys)

% A = sys.a;
% B = sys.b;
% C = sys.c;
% D = sys.d;

freqresp(tf(sys),0)


ClosedLoopModel =linmod('pidCascade_struc2_Analysis|');
Ac = ClosedLoopModel.a;
Bc = ClosedLoopModel.b;
Cc = ClosedLoopModel.c;
Dc = ClosedLoopModel.d;

inputName = [];%cell(length(ClosedLoopModel.InputName),1);
outputName = [];cell(length(ClosedLoopModel.OutputName),1);

for i = 1:length(ClosedLoopModel.InputName)
    
    myString = strsplit(ClosedLoopModel.InputName{i},'/');
    inputName = [inputName;myString(2)];
end


for i = 1:length(ClosedLoopModel.OutputName)
    
    myString = strsplit(ClosedLoopModel.OutputName{i},'/');
    outputName = [outputName;myString(2)];
end

sysc = ss(Ac,Bc,Cc,Dc,'InputName',inputName,'OutputName',outputName);

pc = pole(sysc);

sort(pc)

S11 = sysc('yT1','dT1');
S22 = sysc('yT2','dT2');
S33 = sysc('yTOP','dTOP');

T11 = sysc('yT1','rT1');
T22 = sysc('yT2','rT2');
T33 = sysc('yTOP','rTOP');

sim('pidCascade_struc2_Tuner')

t = t/3600;
%%
figure(1)
clf
subplot(3,1,1)
plot(t,PCT1y,t,PCT1r)
subplot(3,1,2)
plot(t,PCT2y,t,PCT2r)
subplot(3,1,3)
plot(t,PCTOPy,t,PCTOPr)

figure(2)
clf
subplot(3,1,1)
plot(t,PCA1y,t,PCA1r)
subplot(3,1,2)
plot(t,PCA2y,t,PCA2r)
subplot(3,1,3)
plot(t,PCINLy,t,PCINLr)

figure(3)
clf
subplot(3,1,1)
plot(t,z1)
subplot(3,1,2)
plot(t,z2)
subplot(3,1,3)
plot(t,z3)
%%
omega = logspace(-5,1,1000);

absS11 = reshape(abs(freqresp(S11,omega)),1,1000);
absS22 = reshape(abs(freqresp(S22,omega)),1,1000);
absS33 = reshape(abs(freqresp(S33,omega)),1,1000);

absT11 = reshape(abs(freqresp(T11,omega)),1,1000);
absT22 = reshape(abs(freqresp(T22,omega)),1,1000);
absT33 = reshape(abs(freqresp(T33,omega)),1,1000);

figure(4)
clf
subplot(2,1,1)
semilogx(omega,absS11,omega,absT11)
subplot(2,1,2)
semilogx(omega,absS22,omega,absT22)
disp('----------------------')
disp(['Peak S11 = ' num2str(max(absS11))])
disp(['Peak T11 = ' num2str(max(absT11))])
disp('----------------------')
disp(['Peak S22 = ' num2str(max(absS22))])
disp(['Peak T22 = ' num2str(max(absT22))])
disp('----------------------')
disp(['Peak S33 = ' num2str(max(absS33))])
disp(['Peak T33 = ' num2str(max(absT33))])
disp('----------------------')
%%

col1 = [8, 118, 184]/255;
col2 = [36, 162, 14]/255;
col3 = [233, 11, 0]/255;


figure(5)
clf
rect = [0, 0, 12, 15];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[12 15],'PaperPosition',rect)
subplot(3,1,1)
p1 = semilogx(omega,absS11,'-','Color',col1,'LineWidth',1);
hold on
p2 = semilogx(omega,absT11,'--','Color',col3,'LineWidth',1);
L1 = legend([p1 p2],'|\itS\rm|','|\itT\rm|','Location','northwest','Orientation','Vertical','Interpreter','latex');
% legend('boxoff')
% set(L1,'Position',[0.1,0.5+0.12,0.5,0.05])
%text(8,35,'open-loop stable','FontWeight','Bold','FontSize',16,'Color',[1 1 1])  
%text(7.5,29,'open-loop unstable','FontWeight','Bold','FontSize',16,'Color',[0.8 0.8 0.8])  
title('Sensitivity transfer functions for CS2 applied on well A','Interpreter','latex')
xlabel('$\omega$ [rad/s]','Interpreter','latex')
ylabel('Sensitivity','Interpreter','latex')
yticks([0 0.5 1 1.4])
xlim([-5 1])
ylim([0 1.4])

%%ylim([50 51])
subplot(3,1,2)
p3 = semilogx(omega,absS22,'-','Color',col1,'LineWidth',1);
hold on
p4 = semilogx(omega,absT22,'--','Color',col3,'LineWidth',1);
L2 = legend([p3 p4],'|\itS\rm|','|\itT\rm|','Location','northwest','Orientation','Vertical','Interpreter','latex');
uistack(p3,'top')
uistack(p4,'bottom')
yticks([0 0.5 1 1.4])
xlim([-5 1])
ylim([0 1.4])
xlabel('$\omega$ [rad/s]','Interpreter','latex')
ylabel('Sensitivity','Interpreter','latex')
title('Sensitivity transfer functions for CS2 applied on well B','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';

subplot(3,1,3)
p5 = semilogx(omega,absS33,'-','Color',col1,'LineWidth',1);
hold on
p6 = semilogx(omega,absT33,'--','Color',col3,'LineWidth',1);
L3 = legend([p5 p6],'|\itS\rm|','|\itT\rm|','Location','northwest','Orientation','Vertical','Interpreter','latex');
uistack(p5,'top')
uistack(p6,'bottom')
yticks([0 0.5 1 1.4])
xlim([-5 1])
ylim([0 1.4])
xlabel('$\omega$ [rad/s]','Interpreter','latex')
ylabel('Sensitivity','Interpreter','latex')
title('Sensitivity transfer functions for CS4 applied on pipeline-riser','Interpreter','latex')
axs = gca;
axs.TickLabelInterpreter = 'latex';

print -depsc C:\Users\Admin\Dropbox\OptimumSeeking\Manuscript2_revised\Figures\Sensitivities
