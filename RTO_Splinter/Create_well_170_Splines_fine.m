clc
clear all
warning off
close all

load_libSplinter()
%load 'Well_170_Data.mat';
load 'C:\Users\Admin\Dropbox\Journal_Papers\Plantwide\olgaRead\Well_170_Olga2014_fine_Data.mat';

NP = length(Pwh);
NW = length(WGin);

xs = zeros(NP*NW, 2);
P_whds = zeros(NP*NW, 1);
T_whds = zeros(NP*NW, 1);
WG_whds = zeros(NP*NW, 1);
WL_whds = zeros(NP*NW, 1);
k = 1;
for i = 1:NP
   for j = 1:NW
      xs(k,:) = [1e5*Pwh(i) WGin(j)];
      P_whds(k) = P_whd(i,j);
      T_whds(k) = T_whd(i,j);
      WG_whds(k) = WG_whd(i,j);
      WL_whds(k) = WL_whd(i,j);
      k = k + 1;
   end
end


figure(2)
clf
surf(WGin,1e-5*Pwh,WL_whd)
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
zlabel('Oil production rate [kg/s]')
%%
rect = [0, 0, 12, 12];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
figure(3)
clf
surf(WGin,1e-5*Pwh,QL_whd)
axis tight
view(-50,30)
%camlight left
zlabel('Produced oil [STB/day]')
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')

%print -depsc OilFlowRateSpline
%%
figure(4)
clf
surf(WGin,1e-5*Pwh,T_whd,'FaceAlpha',1)
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
zlabel('Temperaure [C]')


%%
figure(6)
clf
surf(WGin,1e-5*Pwh,P_whd,'FaceAlpha',1)
zlabel('P_{wh}')
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
%%

WL_w170 = BSplineBuilder(xs, WL_whds, 4).build();
WL_w170.save('WL_w170_Olga2014_fine.bspline');

WG_w170 = BSplineBuilder(xs, WG_whds, 4).build();
WG_w170.save('WG_w170_Olga2014_fine.bspline');

PWH_w170 = BSplineBuilder(xs, P_whds, 4).build();
PWH_w170.save('PWH_w170_Olga2014_fine.bspline');

TWH_w170 = BSplineBuilder(xs, T_whds, 4).build();
TWH_w170.save('TWH_w170_Olga2014_fine.bspline');

% save('IPR_170_2D_Splines_mass.mat','WL_w170','QL_w170','WG_w170','QG_w170','PWH_w170','TWH_w170')

WL = WL_w170.eval([28.5 1.35])
T1 = TWH_w170.eval([28.5 1.35])
%%
X2 = linspace(min(1e5*Pwh),max(1e5*Pwh),200);
X3 = linspace(min(WGin),max(WGin),200);
%Z = reshape(fnval( QL_w170, {X2, X3}),200,200);
Z = zeros(200,200);
for i=1:200
    for j=1:200
        Z(i,j) = WL_w170.eval([X2(i) X3(j)]);
    end
end
figure(7)
surf(X3,X2,Z,'LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
%daspect([5 5 1])
axis tight
view(-50,30)
camlight left
zlabel('Produced oil [kg/day]')
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
