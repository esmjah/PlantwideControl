clc
clear all
warning off
close all

load_libSplinter()

[FCP_G,FCP_L] = Build_PVT_functions();

load C:\olgaRTO2014\olgaRead\PipelineData_new.mat
%load PipelineData.mat
%load Pipeline_Olga2014_Data.mat

%Pout = 1e-5*Pout;
%PINL = 1e-5*PINL;


NT = length(Tin);
NF = length(Fin);
NW = length(Win);
NP = length(Pout);

Hin = zeros(NT,NF,NW,NP);
Ein = zeros(NT,NF,NW,NP);
T1 = zeros(NT,NF,NW,NP);
T2 = zeros(NT,NF,NW,NP);


PINLs = zeros(NT*NF*NW*NP,1);
TINLs = zeros(NT*NF*NW*NP,1);
TTOPs = zeros(NT*NF*NW*NP,1);
xs = zeros(NT*NF*NW*NP,4);
n = 1;
tic
for m=1:NT
    for i=1:NF
        for j=1:NW
            for k =1:NP
                
                xs(n,:) = [TINL_av(m) Fin(i) Win(j) Pout(k)];
                
                T1(m,i,j,k) = TINL(m,i,j,k); %Tin(m);
                T2(m,i,j,k) = TTOP(m,i,j,k);
                Fin1 = Fin(i);
                W1 = Win(j);
                P1 = PINL(m,i,j,k);
                P2 = Pout(k);
                
                WG = Fin1*W1;
                WO = (1-Fin1)*W1;
                
%                 Hg = fnval(FH_G,{P1,T1(m,i,j,k)}) - fnval(FH_G,{1e5,0});
%                 Ho = fnval(FH_L,{P1,T1(m,i,j,k)}) - fnval(FH_G,{1e5,0});
%                 
%                 Hin(m,i,j,k) = Hg*WG + Ho*WO;
                
                % *************************************
                Cpg = FCP_G.eval([P1 T1(m,i,j,k)]);
                Cpo = FCP_L.eval([P1 T1(m,i,j,k)]);
                Ein(m,i,j,k) = T1(m,i,j,k)*(Cpg*WG + Cpo*WO);
                
                PINLs(n) = PINL(m,i,j,k);
                TINLs(n) = TINL(m,i,j,k);
                TTOPs(n) = TTOP(m,i,j,k);
                
                n = n + 1;
            end
            
        end
    end
end

disp('finished reading values.')

FP_INL = BSplineBuilder(xs,PINLs, 2).build();
FP_INL.save('FP_INL_Olga2014_new.bspline');
disp('finished building FP_INL')

FT_INL = BSplineBuilder(xs,TINLs, 2).build();
FT_INL.save('FT_INL_Olga2014_new.bspline');
disp('finished building FT_INL')

FT_TOP = BSplineBuilder(xs,TTOPs, 2).build();
FT_TOP.save('FT_TOP_Olga2014_new.bspline');
disp('finished building FT_TOP')

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