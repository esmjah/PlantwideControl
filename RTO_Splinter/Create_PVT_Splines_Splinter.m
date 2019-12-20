clc
clear all
close all

load PVT.mat
load_libSplinter()

ZCG = reshape(CP_G',10,10)';
ZCL = reshape(CP_L',10,10)';

ZHG = reshape(H_G',10,10)';
ZHL = reshape(H_L',10,10)';

ZROG = reshape(ROG',10,10)';
ZROL = reshape(ROL',10,10)';

ZRS = reshape(RS',10,10)';


T = linspace(0,200,10)';
P = linspace(1,200,10)';

NP = length(P);
NT = length(T);

xs = zeros(NP*NT,2);
dCG = zeros(NP*NT,1);
dCL = zeros(NP*NT,1);
k = 1;
for i = 1:NP
   for j=1:NT
       xs(k,:) = [P(i) T(j)];
      dCG(k) = ZCG(i,j);
      dCL(k) = ZCL(i,j);
      k = k + 1;
   end
end

FCP_G = BSplineBuilder(xs, dCG, 4).build();

FCP_L = BSplineBuilder(xs, dCL, 4).build();

FCP_G.save('FCP_G_4.bspline');
FCP_L.save('FCP_L_4.bspline');

X = linspace(1,200,100)';
Y = linspace(0,200,100)';

S_CG = zeros(length(X),length(Y));
S_CL = zeros(length(X),length(Y));

for i = 1:length(X)
   for j=1:length(Y)
      S_CG(i,j) = FCP_G.eval([X(i) Y(j)]);
      S_CL(i,j) = FCP_L.eval([X(i) Y(j)]);
   end
end

%%
figure(1)
rect = [0, 0, 18, 12];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[18 12],'PaperPosition',rect)
clf
surf(P,T,ZCG','FaceAlpha',0);
hold on
surf(X, Y, S_CG','LineStyle','none','FaceAlpha',1);
%shading faceted
colormap('gray')
xlim([1 200])
xlabel('PRESSURE [Bara]')
ylim([0 200])
ylabel('TEMPERATURE [C]')
zlabel('GAS SPECIFIC HEAT [J/kg-C]')
view(-45,30)
%%
figure(2)
rect = [0, 0, 18, 12];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[18 12],'PaperPosition',rect)
clf
surf(P,T,ZCL','FaceAlpha',0);
hold on
surf(X, Y, S_CL','LineStyle','none','FaceAlpha',1);
%shading faceted
colormap('gray')
xlim([1 200])
xlabel('PRESSURE [Bara]')
ylim([0 200])
ylabel('TEMPERATURE [C]')
zlabel('OIL SPECIFIC HEAT [J/kg-C]')
view(-45,30)
