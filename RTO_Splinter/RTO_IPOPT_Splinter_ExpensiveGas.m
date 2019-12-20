clc
clear all
close all

load_libSplinter()

%load ('Hgrid.mat');

FTH_TOP = BSpline('FTH_TOP.bspline');
FT_TOP  = BSpline('FT_TOP.bspline');
FP_INL  = BSpline('FP_INL.bspline');
FT_INL  = BSpline('FT_INL.bspline');
FCP_G   = BSpline('FCP_G.bspline');
FCP_L   = BSpline('FCP_L.bspline');

WL_w160   = BSpline('WL_w160.bspline');
WG_w160   = BSpline('WG_w160.bspline');
PWH_w160  = BSpline('PWH_w160.bspline');
TWH_w160  = BSpline('TWH_w160.bspline');

WL_w170   = BSpline('WL_w170.bspline');
WG_w170   = BSpline('WG_w170.bspline');
PWH_w170  = BSpline('PWH_w170.bspline');
TWH_w170  = BSpline('TWH_w170.bspline');

w1 = [];
w2 = [];

global P_sep DeltaP1_manifold DeltaP2_manifold
P_sep = 5.1e5;
DeltaP1_manifold = 5000;%0.6e5;
DeltaP2_manifold = 5000;%0.83e5;

w1.FWG = WG_w160;
w2.FWG = WG_w170;

w1.FWL = WL_w160;
w2.FWL = WL_w170;

w1.FPWH = PWH_w160;
w2.FPWH = PWH_w170;

w1.FTWH = TWH_w160;
w2.FTWH = TWH_w170;


Lambda(1) = 1;
Lambda(2) = 2;
Winj1 = 1.6;
Winj2 = 1.6;
Delta_Pwh1 = 5;
Delta_Pwh2= 5;
Delta_Ptop = 5;
Pin = 23;

  x0         = [Winj1 , Winj2, Delta_Pwh1, Delta_Pwh2, Delta_Ptop, Pin]';  % The starting point.
  options.lb = [0.5, 0.5, 2.2,  2.2, 0.75, 15]';   % Lower bound on the variables.
  options.ub = [2.0, 2.0, 20, 20, 20, 30]';  % Upper bound on the variables.
  options.cl = [0 0];   % Lower bounds on the constraint functions.
  options.cu = [0 inf]; % Upper bounds on the constraint functions.

  % Set the IPOPT options.
  options.ipopt.jac_c_constant        = 'yes';
  options.ipopt.hessian_approximation = 'limited-memory';
  options.ipopt.print_user_options    = 'yes';
  options.ipopt.mu_strategy           = 'adaptive';
  options.ipopt.tol                   = 1e-8;
  options.ipopt.constr_viol_tol       = 1e-8;
  options.ipopt.compl_inf_tol         = 1e-8;
  options.ipopt.dual_inf_tol          = 1e-8;
  

  % The callback functions.
  funcs.objective         = @(x)CostFuncScaled_splinter(x,Lambda,w1.FWL,w2.FWL);
  funcs.constraints       = @(x)NonlinConScaled_splinter( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FTH_TOP);
  funcs.gradient          = @(x)gradient_splinter(x,Lambda,w1.FWL,w2.FWL);
  funcs.jacobian          = @(x)jacobian_splinter( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FTH_TOP);
  funcs.jacobianstructure = @()sparse(ones(2,6));%sparse(tril(ones(2,6)));
%  funcs.iterfunc          = @callback;


%   funcs.gradient          = @gradient;
%   funcs.jacobian          = @jacobian;
%   funcs.jacobianstructure = @jacobian;

[xs info] = ipopt(x0,funcs,options);

x(1:2) = xs(1:2);
x(3:6) = 1e5*xs(3:6);

Winj1 = x(1);
Winj2 = x(2);
Pwh1 = x(3) + x(6);
Pwh2 = x(4) + x(6);
Ptop = x(5) + 5.0816e5;
Pin = x(6);

WL1 = w1.FWL.eval([Pwh1 Winj1]);
WL2 = w2.FWL.eval([Pwh2 Winj2]);

Jopt = -Lambda(1)*(WL1+WL2) + Lambda(2)*(Winj1+Winj2);

disp('***** Optimal solution *****')
disp(['Obj. = ' num2str(Jopt)]);
disp(['Winj1 = ' num2str(x(1))]);
disp(['Winj2 = ' num2str(x(2))]);
disp(['Delta_Pwh1 = ' num2str(x(3))]);
disp(['Delta_Pwh2 = ' num2str(x(4))]);
disp(['Delta_Ptop = ' num2str(x(5))]);
disp(['Pin = ' num2str(Pin)]);
disp(['WL1 = ' num2str(WL1)]);
disp(['WL2 = ' num2str(WL2)]);

[WLtot,ALPHA1,ALPHA2,ALPHA,PT1,PT2,Pin,Ptop,T1,T2,Tin,T_top] = CalculateSetpoints_splinter( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FTH_TOP);

disp('***** Optimal setpoints *****')
disp(['WL = ' num2str(WLtot) ' kg/sec'])

disp(['ALPHA1 = ' num2str(ALPHA1)])
disp(['ALPHA2 = ' num2str(ALPHA2)])
disp(['ALPHA = ' num2str(ALPHA)])
disp(['Pin = ' num2str(Pin)])
disp(['Ptop = ' num2str(Ptop)])
disp(['T1 = ' num2str(T1)])
disp(['T2 = ' num2str(T2)])
disp(['Tin = ' num2str(Tin)])
disp(['Ttop = ' num2str(T_top)])
disp(['PWD1 = ' num2str(PT1)])
disp(['PWD2 = ' num2str(PT2)])

halt

lb = options.lb;
ub = options.ub;
lb(3:6) = 1e5*lb(3:6); 
ub(3:6) = 1e5*ub(3:6); 
%%
X1 = linspace(0.2,3,200);
X3 = 1e5*linspace(2,10,200);
Z = zeros(200,200);
for i=1:200
    for j=1:200
        Z(i,j) = WL_w160.eval([Pin+X3(i) X1(j)]);
    end
end
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
Z = zeros(200,200);
for i=1:200
    for j=1:200
        Z(i,j) = WL_w170.eval([Pin+X4(i) X2(i)]);
    end
end

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

Z = zeros(200,200);
for i=1:200
    for j=1:200
        Z(i,j) = FP_INL.eval([Tin ALPHA X7(i) X5(j)+5e5]);
    end
end

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
        WL1 = w1.FWL.eval([Pin+X3(i) X1(j)]);
        Z1(i,j) = PriceL*(WL1) - PriceG*(X1(j));
    end
end

C1 = zeros(1,length(X1));

for j=1:length(X1)
    WL1 = w1.FWL.eval([Pin+lb(3) X1(j)]);
    C1(j) = PriceL*(WL1) - PriceG*(X1(j));
end

WL1_opt = w1.FWL.eval([Pin+x(3),x(1)]);
J1_opt = PriceL*WL1_opt - PriceG*x(1);
%%
%Z = reshape(fnval( w1.FWL, {Pin, X2, X3}),200,200);
figure(4)
clf
surf(X3,X1,Z1','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
hold on
plot3(lb(3)*ones(size(X3)),X1,C1,'--r','LineWidth',3)
plot3(x(3),x(1),J1_opt,'*k','MarkerSize',10) %/125730
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
        WL2 = w2.FWL.eval([x(6)+X4(i) X2(j)]);
        Z2(i,j) = PriceL*(WL2) - PriceG*(X2(j));
    end
end

C2 = zeros(1,length(X2));

for j=1:length(X2)
    WL2 = w2.FWL.eval([Pin+lb(4) X2(j)]);
    C2(j) = PriceL*(WL2) - PriceG*(X2(j));
end

WL2_opt = w2.FWL.eval([Pin+x(4),x(2)]);
J2_opt = PriceL*WL2_opt - PriceG*x(2);
%%
%Z = reshape(fnval( w1.FWL, {Pin, X2, X3}),200,200);
figure(5)
clf
surf(X4,X2,Z2','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
hold on
plot3(lb(4)*ones(size(X4)),X2,C2,'--r','LineWidth',3)
plot3(x(4),x(2),J2_opt,'*b','MarkerSize',10)
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

rect = [0, 0, 18, 12];
grad = linspace(0.2,0.9,128)';
cmap = [grad, grad, grad];
figure(6)
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[18 12],'PaperPosition',rect)
clf
h = surf(1e-5*X3X4,X1X2,J','FaceColor', 'none','LineStyle','none');
colormap(cmap);
set(h,'FaceLighting','phong',...
      'FaceColor','interp',...
      'BackFaceLighting','lit',...
      'AmbientStrength',0.7,...
      'FaceAlpha',0.8)
light('Position',[-0.002 0 0.01],'Style','infinite')
material dull
hold on
plot3(1e-5*(lb(3)+lb(4))*ones(size(X3X4)),X1X2,C1+C2,'--k','LineWidth',1)
hold on
plot3(1e-5*(x(3)+x(4)),x(1)+x(2),J_opt,'hk','MarkerSize',10,'MarkerFaceColor','k')
%daspect([5 5 1])
axis tight
view(45,30)
%camlight left
zlabel('Obj. function','FontSize',14)
ylabel('u_1+u_2 [kg/s]','FontSize',14)
xlabel('u_3+u_4 [bar]','FontSize',14)

%fileName = 'Obj_Gas_new';
%str = ['print -dpng -r600 ' fileName];
%eval(str);