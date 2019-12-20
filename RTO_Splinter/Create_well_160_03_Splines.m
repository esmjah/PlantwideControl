clc
clear all
warning off
close all

load_libSplinter()
%load 'Well_160_Data.mat';
load 'Well_160_03_Data.mat';

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
      xs(k,:) = [Pwh(i) WGin(j)];
      P_whds(k) = P_whd(i,j);
      T_whds(k) = T_whd(i,j);
      WG_whds(k) = WG_whd(i,j);
      WL_whds(k) = WL_whd(i,j);
      k = k + 1;
   end
end


figure(2)
clf
surf(WGin,Pwh,WL_whd)
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
zlabel('Oil production rate [kg/s]')
%%
rect = [0, 0, 12, 12];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
figure(3)
clf
surf(WGin,Pwh,QL_whd)
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
surf(WGin,Pwh,T_whd,'FaceAlpha',1)
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
zlabel('Temperaure [C]')


%%
figure(6)
clf
surf(WGin,Pwh,P_whd,'FaceAlpha',1)
zlabel('P_{wh}')
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
%%

WL_w160 = BSplineBuilder(xs, WL_whds, 4).build();
WL_w160.save('WL_w160_03_Olga2014_4.bspline');

WG_w160 = BSplineBuilder(xs, WG_whds, 4).build();
WG_w160.save('WG_w160_03_Olga2014_4.bspline');

PWH_w160 = BSplineBuilder(xs, P_whds, 4).build();
PWH_w160.save('PWH_w160_03_Olga2014_4.bspline');

TWH_w160 = BSplineBuilder(xs, T_whds, 4).build();
TWH_w160.save('TWH_w160_03_Olga2014_4.bspline');

% save('IPR_160_2D_Splines_mass.mat','WL_w160','QL_w160','WG_w160','QG_w160','PWH_w160','TWH_w160')

WL = WL_w160.eval([28.5 1.35])
T1 = TWH_w160.eval([28.5 1.35])
%%
X2 = linspace(min(Pwh),max(Pwh),200);
X3 = linspace(min(WGin),max(WGin),200);
%Z = reshape(fnval( QL_w160, {X2, X3}),200,200);
Z = zeros(200,200);
for i=1:200
    for j=1:200
        Z(i,j) = WL_w160.eval([X2(i) X3(j)]);
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

WG_whd_spline = zeros(NP,NW);
WL_whd_spline = zeros(NP,NW);

WG_whd_error = zeros(NP,NW);
WL_whd_error = zeros(NP,NW);

for i = 1:NP
   for j = 1:NW
      WG_whd_spline(i,j) = WG_w160.eval([Pwh(i) WGin(j)]);
      WL_whd_spline(i,j) = WL_w160.eval([Pwh(i) WGin(j)]);
      WG_whd_error(i,j) = WG_whd(i,j) - WG_whd_spline(i,j);
      WL_whd_error(i,j) = WL_whd(i,j) - WL_whd_spline(i,j);
   end
end


figure(8)
clf
surf(WGin,Pwh,WL_whd_error)
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
zlabel('Oil production rate error [kg/s]')


figure(9)
clf
surf(WGin,Pwh,WG_whd_error)
xlabel('Gas injection rate [kg/s]')
ylabel('Well-head pressure [bar]')
zlabel('Gas production rate error [kg/s]')