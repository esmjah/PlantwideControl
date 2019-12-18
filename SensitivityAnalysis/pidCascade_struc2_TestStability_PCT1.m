clc
clear all

format long g

load jacobians_struc2_z5.mat

close all
warning off

inputNames = {'inj1','inj2','zA1','zA2','zINL'}';
outputNames = {'yA1','yT1','yA2','yT2','yINL','yTOP'}';

sys = ss(A,B,C,D,'InputName',inputNames,'OutputName',outputNames);

s = tf('s');
%sys = sys_in;% = minreal(sys_in);

%pole(sys)

%freqresp(tf(sys),0)

load('TuningValues.mat')

% Tubing Well #1
PCT1Kc = 0.016/40e5;
PCT1Ti = 100;
PCT1Tf = 0;
PCT1 = PCT1Kc*(1+1/(s*PCT1Ti))*(40e5/(1+s*PCT1Tf));

save('TuningValues.mat','PCINLKc','PCINLTi','PCINLTf','PCTOPKc','PCTOPTi','PCTOPTf','PCA1Kc','PCA1Ti','PCA1Tf','PCT1Kc','PCT1Ti','PCT1Tf','PCA2Kc','PCA2Ti','PCA2Tf','PCT2Kc','PCT2Ti','PCT2Tf')

ClosedLoopModel =linmod('pidCascade_struc2_Analysis_PCA1');
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

%pc = pole(sysc);

%sort(pc)

GA1 = sysc('yA1','z1');

ClosedLoopModel =linmod('pidCascade_struc2_Analysis_PCT1');
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

PCA1 = PCA1Kc*(1+1/(s*PCA1Ti))*(1/(1+s*PCA1Tf));

LA1 = series(PCA1,GA1);
GCA1 = feedback(LA1,1);
KSA1 = feedback(PCA1,GA1);
mar_A1 = allmargin(LA1);
disp('PCA1 margins:')
disp(mar_A1)

GT1 = sysc('yT1','rA1');
LT1 = series(PCT1,GT1);
GCT1 = feedback(LT1,1);

mar_T1 = allmargin(LT1);
disp('PCT1 margins:')
disp(mar_T1)

figure(1)
margin(LA1)
figure(2)
margin(LT1)
