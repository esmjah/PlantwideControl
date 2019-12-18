clc
clear all

format long g
load jacobians_struc2_zOpt.mat

warning off

inputNames = {'inj1','inj2','zA1','zA2','zINL'}';
outputNames = {'yA1','yT1','yA2','yT2','yINL','yTOP'}';

sys = ss(A,B,C,D,'InputName',inputNames,'OutputName',outputNames);

s = tf('s');
%sys = sys_in;% = minreal(sys_in);

%pole(sys)
%freqresp(tf(sys),0)

load('TuningValues.mat')


PCA2Kc = 3*(-1e-5);
PCA2Ti = 600;
PCA2Tf = 0;
PCA2 = PCA2Kc*(1+1/(s*PCA2Ti))*(1/(1+s*PCA2Tf));

save('TuningValues.mat','PCINLKc','PCINLTi','PCINLTf','PCTOPKc','PCTOPTi','PCTOPTf','PCA1Kc','PCA1Ti','PCA1Tf','PCT1Kc','PCT1Ti','PCT1Tf','PCA2Kc','PCA2Ti','PCA2Tf','PCT2Kc','PCT2Ti','PCT2Tf')


ClosedLoopModel =linmod('pidCascade_struc2_Analysis_PCA2');
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

GA2 = sysc('yA2','z2');

LA2 = series(PCA2,GA2);
GCA2 = feedback(LA2,1);
KSA2 = feedback(PCA2,GA2);
mar_A2 = allmargin(LA2);
disp('PCA2 margins:')
disp(mar_A2)