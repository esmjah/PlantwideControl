clc
clear all
close all

load Pipeline_Splines.mat
load PVT_Splines.mat
[w1,w2] = LoadWellData_2D_mass(160,170);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.5, 0.5, 3e5,  3e5, 1.5e5, 20e5]';
ub = [2.0, 2.0, 20e5, 20e5, 20e5, 30e5]';

options = optimoptions('fmincon','Algorithm','interior-point','MaxIter',300,'TolFun',1e-20,'TolCon',1e-20,'TolX',1e-20,'MaxFunEvals',1e6,'Display','iter');


Lambda(1) = 1;
Lambda(2) = 0.25;
f = @(x)CostFunc_2D_mass(x,Lambda,w1.FWL,w2.FWL);
nonlcon = @(x)NonlinCon_2D_mass( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FH_G,FH_L,FTH_TOP);

Winj1 = 1;
Winj2 = 1;
Delta_Pwh1 = 5e5;
Delta_Pwh2= 5e5;
Delta_Ptop = 5e5;
Pin = 23e5;
x0 = [Winj1 , Winj2, Delta_Pwh1, Delta_Pwh2, Delta_Ptop, Pin]';

[x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

Winj1 = x(1);
Winj2 = x(2);
Pwh1 = x(3) + x(6);
Pwh2 = x(4) + x(6);
Ptop = x(5) + 5e5;
Pin = x(6);

WL1 = fnval(w1.FWL,{Pwh1,Winj1});
WL2 = fnval(w2.FWL,{Pwh2,Winj2});

disp('***** Optimal solution *****')
disp(['Obj. = ' num2str(fval)]);
disp(['Winj1 = ' num2str(x(1))]);
disp(['Winj2 = ' num2str(x(2))]);
disp(['Delta_Pwh1 = ' num2str(x(3))]);
disp(['Delta_Pwh2 = ' num2str(x(4))]);
disp(['Delta_Ptop = ' num2str(x(5))]);
disp(['Pin = ' num2str(Pin)]);
disp(['WL1 = ' num2str(WL1)]);
disp(['WL2 = ' num2str(WL2)]);

[WLtot,ALPHA1,ALPHA2,ALPHA,PT1,PT2,Pin,T1,T2,Tin,T_top] = CalculateSetpoints_2D_mass( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FH_G,FH_L,FTH_TOP);

disp('***** Optimal setpoints *****')
disp(['WL = ' num2str(WLtot) ' kg/sec'])

disp(['ALPHA1 = ' num2str(ALPHA1)])
disp(['ALPHA2 = ' num2str(ALPHA2)])
disp(['ALPHA = ' num2str(ALPHA)])
disp(['Pin = ' num2str(Pin)])
disp(['T1 = ' num2str(T1)])
disp(['T2 = ' num2str(T2)])
disp(['Tin = ' num2str(Tin)])
disp(['Ttop = ' num2str(T_top)])
disp(['PWD1 = ' num2str(PT1)])
disp(['PWD2 = ' num2str(PT2)])

%%
X1 = linspace(0.2,3,200);
X3 = 1e5*linspace(2,10,200);
Z = reshape(fnval( w1.FWL, {Pin+X3, X1}),200,200);
figure(1)
clf
surf(X3,X1,Z','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
%daspect([5 5 1])
axis tight
view(-50,30)
camlight left
zlabel('Produced oil [kg/sec]')
ylabel('Injected gas [kg/sec]')
xlabel('Valve pressure drop [Pa]')
%%
X2 = linspace(0.2,3,200);
X4 = 1e5*linspace(2,10,200);
Z = reshape(fnval( w2.FWL, {Pin+X4, X2}),200,200);
figure(2)
clf
surf(X4,X2,Z','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
%daspect([5 5 1])
axis tight
view(-50,30)
camlight left
zlabel('Produced oil [kg/sec]')
ylabel('Injected gas [kg/sec]')
xlabel('Valve pressure drop [Pa]')

%%
X7 = linspace(22,40,200);
X5 = 1e5*linspace(0.1,9,200);
Z = reshape(fnval(FP_INL,{Tin,ALPHA,X7,X5+5e5}),200,200);
figure(3)
clf
surf(X7,X5,Z','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
%daspect([5 5 1])
axis tight
view(-120,30)
%camlight right
zlabel('Manifold pressure [Pa]')
xlabel('Mass flow [kg/sec]')
ylabel('Vlave pressure drop [Pa]')
%%
PriceL = Lambda(1);
PriceG = Lambda(2);
X1 = linspace(0.2,3,200);
X3 = 1e5*linspace(1,10,200);
Z1 = zeros(length(X3),length(X1));
for i=1:length(X3)
    for j=1:length(X1)
        %XX = [X3(j),X3(j),X2(i),X2(i),x(5),x(6)]';
        %Z(i,j) = feval(f,XX);
        WL1 = fnval(w1.FWL,{Pin+X3(i),X1(j)});
        Z1(i,j) = -PriceL*(WL1) + PriceG*(X1(j));
    end
end

C1 = zeros(1,length(X1));

for j=1:length(X1)
    WL1 = fnval(w1.FWL,{Pin+lb(3),X1(j)});
    C1(j) = -PriceL*(WL1) + PriceG*(X1(j));
end

WL1_opt = fnval(w1.FWL,{Pin+x(3),x(1)});
J1_opt = -PriceL*WL1_opt + PriceG*x(1);
%%
%Z = reshape(fnval( w1.FWL, {Pin, X2, X3}),200,200);
figure(4)
clf
surf(X3,X1,Z1','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
hold on
plot3(lb(3)*ones(size(X3)),X1,C1,'--r','LineWidth',3)
plot3(x(3),x(1),J1_opt,'*y','MarkerSize',10) %/125730
%daspect([5 5 1])
axis tight
view(260,30)
%camlight left
zlabel('J [-1*USD]')
ylabel('Injected gas [kg/sec]')
xlabel('Valve pressure drop [Pa]')
%%
PriceL = Lambda(1);
PriceG = Lambda(2);
X2 = linspace(0.2,3,200);
X4 = 1e5*linspace(1,10,200);
Z2 = zeros(length(X4),length(X2));
for i=1:length(X4)
    for j=1:length(X2)
        %XX = [X3(j),X3(j),X2(i),X2(i),x(5),x(6)]';
        %Z(i,j) = feval(f,XX);
        WL2 = fnval(w2.FWL,{x(6)+X4(i),X2(j)});
        Z2(i,j) = -PriceL*(WL2) + PriceG*(X2(j));
    end
end

C2 = zeros(1,length(X2));

for j=1:length(X2)
    WL2 = fnval(w2.FWL,{Pin+lb(4),X2(j)});
    C2(j) = -PriceL*(WL2) + PriceG*(X2(j));
end

WL2_opt = fnval(w2.FWL,{Pin+x(4),x(2)});
J2_opt = -PriceL*WL2_opt + PriceG*x(2);
%%
%Z = reshape(fnval( w1.FWL, {Pin, X2, X3}),200,200);
figure(5)
clf
surf(X4,X2,Z2','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
hold on
plot3(lb(4)*ones(size(X4)),X2,C2,'--r','LineWidth',3)
plot3(x(4),x(2),J2_opt,'*y','MarkerSize',10)
%daspect([5 5 1])
axis tight
view(250,30)
%camlight left
zlabel('J [-1*USD]')
ylabel('Injected gas [Sm3/day]')
xlabel('Valve pressure drop [Pa]')
%%
X1X2 = X1+X2;
X3X4 = X3+X4;
J = Z1+Z2;
J_opt = J1_opt + J2_opt;
figure(6)
rect = [0, 0, 18, 12];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[18 12],'PaperPosition',rect)
clf
surf(X3X4,X1X2,J','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
hold on
plot3((lb(3)+lb(4))*ones(size(X3X4)),X1X2,C1+C2,'-r','LineWidth',2)
plot3(x(3)+x(4),x(1)+x(2),J_opt,'*y','MarkerSize',10)
%daspect([5 5 1])
axis tight
view(-120,30)
%camlight left
zlabel('Obj. function')
ylabel('u_1+u_2 [kg/s]')
xlabel('u_3+u_4 [Pa]')

fileName = 'Cost_Cheap_Gas';
str = ['print -dpng -r600 Figures\' fileName];
%eval(str);