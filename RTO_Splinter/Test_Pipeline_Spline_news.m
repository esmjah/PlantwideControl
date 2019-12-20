clc
clear all
warning off
close all

load_libSplinter()

load C:\olgaRTO2014\olgaRead\PipelineData_new.mat

NT = length(Tin);
NF = length(Fin);
NW = length(Win);
NP = length(Pout);

FT_TOP  = BSpline('FT_TOP_Olga2014_new.bspline');
FP_INL  = BSpline('FP_INL_Olga2014_new.bspline');
FT_INL  = BSpline('FT_INL_Olga2014_new.bspline');
FCP_G   = BSpline('FCP_G_4.bspline');
FCP_L   = BSpline('FCP_L_4.bspline');
%%
PINL_test = FP_INL.eval([95.62 7.7 33 5.73])

PINL_jac_test = FP_INL.eval_jacobian([95.62 7.7 33 5.73])
%%
m0 = 4;
i0 = 4;
rect = [0, 0, 15, 15];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
TTOP1 = reshape(TTOP(m0,i0,:,:),NW,NP)';
figure(1)
surf(Pout,Win,TTOP1')
axis tight
view(50,30)
%camlight left
zlabel('TM TOP []')
ylabel('inlet mass flow [kg/sec]')
xlabel('top pressure [bar]')

%print -depsc ManifoldPressureSplines
%%
X3 = linspace(min(Win),max(Win),200);
X4 = linspace(min(Pout),max(Pout),200);
%Z = reshape(fnval( QL_w160, {X2, X3}),200,200);
Z = zeros(200,200);
for n=1:200
    for j=1:200
        Z(n,j) = FP_INL.eval([TINL_av(m0) Fin(i0) X3(j) X4(n)]);
    end
end
%%
figure(2)
clf
rect = [0, 0, 15, 15];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
s1 = surf(X4,X3,Z','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
%daspect([5 5 1])
hold on
P1 = reshape(PINL(m0,i0,:,:),NW,NP)';
s2 = surf(Pout,Win,P1');
alpha(s2,0)
axis tight
view(-50,30)
%camlight left
zlabel('manifold pressure [bar]')
ylabel('total flow rate [kg/s]')
xlabel('top pressure [bar]')
%%
Error = zeros(NP,NW);
for n = 1:NP
    for j = 1:NW
        Error(n,j) = PINL(m0,i0,j,n) - FP_INL.eval([TINL_av(m0) Fin(i0) Win(j) Pout(n)]);
    end
end
figure(3)
clf
rect = [0, 0, 15, 15];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
surf(Pout,Win,Error')
axis tight
view(-50,30)
%camlight left
zlabel('modeling error [bar]')
ylabel('total flow rate [kg/s]')
xlabel('top pressure [bar]')
%%
k0 = 7;
X2 = linspace(min(Fin),max(Fin),100);
X3 = linspace(min(Win),max(Win),200);
%Z = reshape(fnval( QL_w160, {X2, X3}),200,200);
Z = zeros(100,200);
for i=1:100
    for j=1:200
        Z(i,j) = FP_INL.eval([TINL_av(m0) X2(i) X3(j) Pout(k0)]);
    end
end
%%
figure(4)
clf
rect = [0, 0, 15, 15];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
s1 = surf(X2,X3,Z','LineStyle','none','FaceColor','interp','EdgeColor','black','FaceLighting','phong');
%daspect([5 5 1])
hold on
P2 = zeros(NF,NW);
for i = 1:NF
    for j = 1:NW
        P2(i,j) = PINL(m0,i,j,k0);
    end
end
s2 = surf(Fin,Win,P2');
alpha(s2,0)
axis tight
view(-50,30)
%camlight left
zlabel('manifold pressure [bar]')
ylabel('total flow rate [kg/s]')
xlabel('alpha [-]')
%%
Error = zeros(NF,NW);
for i = 1:NF
    for j = 1:NW
        Error(i,j) = PINL(m0,i,j,k0) - FP_INL.eval([TINL_av(m0) Fin(i) Win(j) Pout(k0)]);
    end
end
figure(5)
clf
rect = [0, 0, 15, 15];
set(gcf,'Color',[1,1,1],'PaperUnits','centimeters','PaperSize',[15 15],'PaperPosition',rect)
surf(Fin,Win,Error')
axis tight
view(-50,30)
%camlight left
zlabel('modeling error [bar]')
ylabel('total flow rate [kg/s]')
xlabel('alpha [-]')
%%
