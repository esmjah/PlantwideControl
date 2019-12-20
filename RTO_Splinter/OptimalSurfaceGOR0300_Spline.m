clc
clear all
warning off
close all

load_libSplinter()
%load 'Well_165_Data.mat';
load 'OptimalSurface_GOR0300.mat';

NW1 = length(WGin1);
NW2 = length(WGin2);

xs = zeros(NW1*NW2, 2);
Js = zeros(NW1*NW2, 1);

GORm1s = zeros(NW1*NW2, 1);
GORm2s = zeros(NW1*NW2, 1);
Pwh1s = zeros(NW1*NW2, 1);
Pwh2s = zeros(NW1*NW2, 1);
Pinls = zeros(NW1*NW2, 1);
Zwh1s = zeros(NW1*NW2, 1);
Zwh2s = zeros(NW1*NW2, 1);
Ztops = zeros(NW1*NW2, 1);
WLwh1s = zeros(NW1*NW2, 1);
WLwh2s = zeros(NW1*NW2, 1);
WGwh1s = zeros(NW1*NW2, 1);
WGwh2s = zeros(NW1*NW2, 1);

k = 1;

J_min = min(min(J));
disp(['Mininum cost is ' num2str(J_min)])

for j = 1:NW2
    for i = 1:NW1
        if(J(j,i)==J_min)
            n = (j-1)*NW1 + i;
            WG1_min = WGin1(i);
            WG2_min = WGin2(j);
            disp(['The minimum happens at WGin1_min=' num2str(WG1_min) ' ,WGin2_min=' num2str(WG2_min)]);
            disp(['n=' num2str(n)])
        end
        xs(k,:) = [WGin1(i) WGin2(j)];
        Js(k) = J(j,i);
        GORm1s(k) = GORm1(j,i);
        GORm2s(k) = GORm2(j,i);
        Pwh1s(k) = Pwh1(j,i);
        Pwh2s(k) = Pwh2(j,i);
        Pinls(k) = Pinl(j,i);
        Zwh1s(k) = Zwh1(j,i);
        Zwh2s(k) = Zwh2(j,i);
        Ztops(k) = Ztop(j,i);
        WLwh1s(k) = WLwh1(j,i);
        WLwh2s(k) = WLwh2(j,i);
        WGwh1s(k) = WGwh1(j,i);
        WGwh2s(k) = WGwh2(j,i);
        k = k + 1;
    end
end


figure(1)
clf
surf(WGin1,WGin2,J)
xlabel('Gas injection rate 1 [kg/s]')
ylabel('Gas injection rate 2 [kg/s]')
zlabel('J')
%%
J_GOR0300 = BSplineBuilder(xs, Js, 3).build();
J_GOR0300.save('OptimalSurface_GOR0300.bspline');
GORm1_GOR0300 = BSplineBuilder(xs, GORm1s, 3).build();
GORm2_GOR0300 = BSplineBuilder(xs, GORm2s, 3).build();
Pwh1_GOR0300 = BSplineBuilder(xs, 1e-5*Pwh1s, 3).build();
Pwh2_GOR0300 = BSplineBuilder(xs, 1e-5*Pwh2s, 3).build();
Pinl_GOR0300 = BSplineBuilder(xs, Pinls, 3).build();
Zwh1_GOR0300 = BSplineBuilder(xs, Zwh1s, 3).build();
Zwh2_GOR0300 = BSplineBuilder(xs, Zwh2s, 3).build();
Ztop_GOR0300 = BSplineBuilder(xs, Ztops, 3).build();
WLwh1_GOR0300 = BSplineBuilder(xs, WLwh1s, 3).build();
WLwh2_GOR0300 = BSplineBuilder(xs, WLwh2s, 3).build();
WGwh1_GOR0300 = BSplineBuilder(xs, WGwh1s, 3).build();
WGwh2_GOR0300 = BSplineBuilder(xs, WGwh2s, 3).build();

% save('IPR_165_2D_Splines_mass.mat','WL_w165','QL_w165','WG_w165','QG_w165','PWH_w165','TWH_w165')

J_eval = J_GOR0300.eval([WG1_min WG2_min])

%%
Error = zeros(NW1,NW2);
for j = 1:NW2
    for i = 1:NW1
        Error(i,j) = J(j,i) - J_GOR0300.eval([WGin1(i) WGin2(j)]);
    end
end
figure(2)
clf
rect = [0, 0, 15, 15];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
surf(WGin1,WGin2,Error')
axis tight
view(-50,30)
%camlight left
zlabel('modeling error')
xlabel('Gas injection rate 1 [kg/s]')
ylabel('Gas injection rate 1 [kg/s]')

%%
  x0         = [WG1_min , WG2_min]';  % The starting point.
  options.lb = [0.7, 1.1]';   % Lower bound on the variables.
  options.ub = [0.9, 1.4]';  % Upper bound on the variables.

  % Set the IPOPT options.
  options.ipopt.jac_c_constant        = 'yes';
  options.ipopt.hessian_approximation = 'limited-memory';
  options.ipopt.print_user_options    = 'yes';
  options.ipopt.mu_strategy           = 'adaptive';
  options.ipopt.tol                   = 1e-12;
  options.ipopt.constr_viol_tol       = 1e-12;
  options.ipopt.compl_inf_tol         = 1e-12;
  options.ipopt.dual_inf_tol          = 1e-12;

  % The callback functions.
  funcs.objective         = @(x)objective (x,J_GOR0300);
  funcs.gradient          = @(x)gradient (x,J_GOR0300);

[xs info] = ipopt(x0,funcs,options);

X_opt = xs
J_opt = J_GOR0300.eval([xs(1) xs(2)]);
GORm1_opt = GORm1_GOR0300.eval([xs(1) xs(2)]);
GORm2_opt = GORm2_GOR0300.eval([xs(1) xs(2)]);
Pwh1_opt = Pwh1_GOR0300.eval([xs(1) xs(2)]);
Pwh2_opt = Pwh2_GOR0300.eval([xs(1) xs(2)]);
Pinl_opt = Pinl_GOR0300.eval([xs(1) xs(2)]);
Zwh1_opt = Zwh1_GOR0300.eval([xs(1) xs(2)]);
Zwh2_opt = Zwh2_GOR0300.eval([xs(1) xs(2)]);
Ztop_opt = Ztop_GOR0300.eval([xs(1) xs(2)]);
WLwh1_opt = WLwh1_GOR0300.eval([xs(1) xs(2)]);
WLwh2_opt = WLwh2_GOR0300.eval([xs(1) xs(2)]);
WGwh1_opt = WGwh1_GOR0300.eval([xs(1) xs(2)]);
WGwh2_opt = WGwh2_GOR0300.eval([xs(1) xs(2)]);

disp(['X_opt: ' num2str(X_opt')])
disp(['J_opt: ' num2str(J_opt)])
disp(['GORm1_opt: ' num2str(GORm1_opt)])
disp(['GORm2_opt: ' num2str(GORm2_opt)])
disp(['Pwh1_opt: ' num2str(Pwh1_opt)])
disp(['Pwh2_opt: ' num2str(Pwh2_opt)])
disp(['Pinl_opt: ' num2str(Pinl_opt)])
disp(['Zwh1_opt: ' num2str(Zwh1_opt)])
disp(['Zwh2_opt: ' num2str(Zwh2_opt)])
disp(['Ztop_opt: ' num2str(Ztop_opt)])
disp(['WLwh1_opt: ' num2str(WLwh1_opt)])
disp(['WLwh2_opt: ' num2str(WLwh2_opt)])
disp(['WGwh1_opt: ' num2str(WGwh1_opt)])
disp(['WGwh2_opt: ' num2str(WGwh2_opt)])

disp('*******************************')
J_jac = J_GOR0300.eval_jacobian([xs(1) xs(2)]);
GORm1_jac = GORm1_GOR0300.eval_jacobian([xs(1) xs(2)]);
GORm2_jac = GORm2_GOR0300.eval_jacobian([xs(1) xs(2)]);
Pwh1_jac = Pwh1_GOR0300.eval_jacobian([xs(1) xs(2)]);
Pwh2_jac = Pwh2_GOR0300.eval_jacobian([xs(1) xs(2)]);
Pinl_jac = Pinl_GOR0300.eval_jacobian([xs(1) xs(2)]);
Zwh1_jac = Zwh1_GOR0300.eval_jacobian([xs(1) xs(2)]);
Zwh2_jac = Zwh2_GOR0300.eval_jacobian([xs(1) xs(2)]);
Ztop_jac = Ztop_GOR0300.eval_jacobian([xs(1) xs(2)]);
WLwh1_jac = WLwh1_GOR0300.eval_jacobian([xs(1) xs(2)]);
WLwh2_jac = WLwh2_GOR0300.eval_jacobian([xs(1) xs(2)]);
WGwh1_jac = WGwh1_GOR0300.eval_jacobian([xs(1) xs(2)]);
WGwh2_jac = WGwh2_GOR0300.eval_jacobian([xs(1) xs(2)]);

disp(['J_jac: ' num2str(J_jac)])
disp(['GORm1_jac: ' num2str(GORm1_jac)])
disp(['GORm2_jac: ' num2str(GORm2_jac)])
disp(['Pwh1_jac: ' num2str(Pwh1_jac)])
disp(['Pwh2_jac: ' num2str(Pwh2_jac)])
disp(['Pinl_jac: ' num2str(Pinl_jac)])
disp(['Zwh1_jac: ' num2str(Zwh1_jac)])
disp(['Zwh2_jac: ' num2str(Zwh2_jac)])
disp(['Ztop_jac: ' num2str(Ztop_jac)])
disp(['WLwh1_jac: ' num2str(WLwh1_jac)])
disp(['WLwh2_jac: ' num2str(WLwh2_jac)])
disp(['WGwh1_jac: ' num2str(WGwh1_jac)])
disp(['WGwh2_jac: ' num2str(WGwh2_jac)])
%%
X1 = linspace(min(WGin1),max(WGin1),200);
X2 = linspace(min(WGin2),max(WGin2),200);
%Z = reshape(fnval( QL_w165, {X2, X3}),200,200);
Z = zeros(200,200);
for i=1:200
    for j=1:200
        Z(i,j) = J_GOR0300.eval([X1(i) X2(j)]);
    end
end
%%
figure(3)
clf
s1 = surf(X1,X2,Z','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
hold on
s2 = surf(WGin1,WGin2,J);
alpha(s2,0)
axis tight
view(120,60)
%camlight left
plot3(X_opt(1),X_opt(2),J_opt,'hk','MarkerSize',10,'MarkerFaceColor','g')
xlabel('Gas injection rate 1 [kg/s]')
ylabel('Gas injection rate 2 [kg/s]')
zlabel('J')
%%

function f = objective (x,J_spline)
  WG1 = x(1);
  WG2 = x(2);
  f = J_spline.eval([WG1 WG2]);
end
function g = gradient (x,J_spline)
  WG1 = x(1);
  WG2 = x(2);
  g = J_spline.eval_jacobian([WG1 WG2]);
end