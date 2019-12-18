clc
clear all

format long g
load jacobians.mat

sys = ss(A,B,C,D);

pole(sys)

freqresp(tf(sys),0)

PCT1Kc = 0.004e-5;
PCT1Ti = 900;
PCT1Tf = 900;


PCA1Kc = -1e-5;
PCA1Ti = 600;
PCA1Tf = 1;

PCT2Kc = 0.004e-5;
PCT2Ti = 900;
PCT2Tf = 900;


PCA2Kc = -1e-5;
PCA2Ti = 600;
PCA2Tf = 1;


PCINLKc = -0.2e-5;
PCINLTi = 900;
PCINLTf = 1;



ClosedLoopModel =linmod('pidCascadeAnalysis');
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

T11 = sysc('yT1','rT1');
T22 = sysc('yT2','rT2');

sim('pidCascadeTuner')

t = t/3600;
%%
figure(1)
clf
subplot(3,1,1)
plot(t,PCT1y,t,PCT1r)
subplot(3,1,2)
plot(t,PCT2y,t,PCT2r)
subplot(3,1,3)
plot(t,PCINLy,t,PCINLr)

figure(2)
clf
subplot(2,1,1)
plot(t,PCA1y,t,PCA1r)
subplot(2,1,2)
plot(t,PCA2y,t,PCA2r)

figure(3)
clf
subplot(3,1,1)
plot(t,z1)
subplot(3,1,2)
plot(t,z2)
subplot(3,1,3)
plot(t,z3)
%%
omega = logspace(-5,-1,1000);

absS11 = reshape(abs(freqresp(S11,omega)),1,1000);
absS22 = reshape(abs(freqresp(S22,omega)),1,1000);

absT11 = reshape(abs(freqresp(T11,omega)),1,1000);
absT22 = reshape(abs(freqresp(T22,omega)),1,1000);

figure(4)
clf
subplot(2,1,1)
semilogx(omega,absS11,omega,absT11)
subplot(2,1,2)
semilogx(omega,absS22,omega,absT22)

disp(['Peak S11 = ' num2str(max(absS11))])
disp(['Peak T11 = ' num2str(max(absT11))])
disp(['Peak S22 = ' num2str(max(absS22))])
disp(['Peak T22 = ' num2str(max(absT22))])
