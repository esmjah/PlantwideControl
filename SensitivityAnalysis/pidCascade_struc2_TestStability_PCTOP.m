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

PCTOPKc = 0.4/20e5;
PCTOPTi = 300;
PCTOPTf = 0;
PCTOP = PCTOPKc*(1+1/(s*PCTOPTi))*(20e5/(1+s*PCTOPTf));
%G_TOP = sys('yTOP','zINL');

save('TuningValues.mat','PCINLKc','PCINLTi','PCINLTf','PCTOPKc','PCTOPTi','PCTOPTf','PCA1Kc','PCA1Ti','PCA1Tf','PCT1Kc','PCT1Ti','PCT1Tf','PCA2Kc','PCA2Ti','PCA2Tf','PCT2Kc','PCT2Ti','PCT2Tf')

ClosedLoopModel =linmod('pidCascade_struc2_Analysis_PCINL');
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


GINL = sysc('yINL','z3');


ClosedLoopModel =linmod('pidCascade_struc2_Analysis_PCTOP');
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

GTOP = sysc('yTOP','rINL');

PCINL = PCINLKc*(1+1/(s*PCINLTi))*(1/(1+s*PCINLTf));
LINL = series(PCINL,GINL);
GCINL = feedback(LINL,1);
KSINL = feedback(PCINL,GINL);
mar_INL = allmargin(LINL);
disp('PCINL margins:')
disp(mar_INL)

LTOP = series(PCTOP,GTOP);
GCTOP = feedback(LTOP,1);
mar_TOP = allmargin(LTOP);
disp('PCTOP margins:')
disp(mar_TOP)

figure(1)
clf
margin(LINL)
figure(2)
clf
margin(LTOP)
