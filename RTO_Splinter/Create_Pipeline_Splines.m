clc
clear all
warning off
close all

load_libSplinter()

[FCP_G,FCP_L] = Build_PVT_functions();

%load C:\olgaRTO2014\olgaRead\PipelineData.mat
%load PipelineData.mat
load Pipeline_Olga2014_Data.mat

Pout = 1e-5*Pout;
PINL = 1e-5*PINL;

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
                
                xs(n,:) = [Tin(m) Fin(i) Win(j) Pout(k)];
                
                T1(m,i,j,k) = Tin(m);
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

FP_INL = BSplineBuilder(xs,PINLs, 4).build();
FP_INL.save('FP_INL_Olga2014_4.bspline');
disp('finished building FP_INL')

FT_INL = BSplineBuilder(xs,TINLs, 4).build();
FT_INL.save('FT_INL_Olga2014_4.bspline');
disp('finished building FT_INL')

FT_TOP = BSplineBuilder(xs,TTOPs, 4).build();
FT_TOP.save('FT_TOP_Olga2014_4.bspline');
disp('finished building FT_TOP')

stophere

% Making T2 as a function of H1, T1, P2
Hin = Ein;
ScalingH = 2e5;
Xin = Hin/ScalingH;
NH = 200;

Hmax = max(max(max(max(Hin))));
Hmin = min(min(min(min(Hin))));

Xmax = max(max(max(max(Xin))));
Xmin = min(min(min(min(Xin))));

Tmin = min(Tin);
Tmax = max(Tin);
H_grid = linspace(Hmin,Hmax,NH);
X_grid = linspace(Xmin,Xmax,NH);
T_grid = linspace(Tmin,Tmax,length(Tin));

iter2 = 1;
x1s = zeros(NT*NW,2);
x2s = zeros(NP*length(T_grid)*NH,3);
dT2_3D = zeros(NP*length(T_grid)*NH,1);
for k=1:NP
    
    dX1p = zeros(NT*NW,1);
    
    T1p = zeros(NT,NW);
    T2p = zeros(NT,NW);
    X1p = zeros(NT,NW);
    iter1 = 1;
    for m=1:NT
        for j=1:NW
            T1p(m,j) = T1(m,4,j,k);
            X1p(m,j) = Xin(m,4,j,k);
            T2p(m,j) = T2(m,4,j,k);
            x1s(iter1,:) = [T1p(m,j) X1p(m,j)];
            dX1p(iter1) = T2p(m,j);
            iter1 = iter1 + 1;
        end
    end
    rbf = RadialBasisFunction(x1s,dX1p, RadialBasisFunctionType.Thin_plate_spline);
    
    clear dX1p
    %[fitresult, gof] = FitT2_H1T1(T1p, H1p, T2p,i);
    T2p_fit = zeros(NT,NH);
    
    for m = 1:length(T_grid)
        for n = 1:NH
            T2_sample = rbf.eval([T_grid(m) X_grid(n)]);
            T2p_fit(m,n) = T2_sample;
            x2s(iter2,:) = [T_grid(m) H_grid(n) Pout(k)];
            dT2_3D(iter2)= T2_sample;  %T2_3D(m,n,i) = T2_sample;
            iter2 = iter2 + 1;
        end
    end
    
    clear rbf
    
    figure( 'Name', ['thin plate spline P=' num2str(Pout(k))]);
    plot3(reshape(T1p,NT*NW,1), reshape(X1p,NT*NW,1),reshape(T2p,NT*NW,1),'o','MarkerFaceColor','b')
    hold on
    surf(T_grid,X_grid,T2p_fit','LineStyle','none')
    grid
    ylabel('H_1 [J/s]')
    xlabel('T_1 [C]')
    zlabel('T_2 [C]')
    
end


disp('Building data tables completed.')

FTH_TOP = BSplineBuilder(x2s,dT2_3D, 4).build();
FTH_TOP.save('FTH_TOP.bspline');

save('Hgrid.mat','H_grid');



%FT_INL.save('FT_INL_Olga2014.bspline');

toc
disp('Bsplines are created.')

disp('Test started.')
tic
Error = zeros(NT,NF,NW,NP);
Hin_est = zeros(NT,NF,NW,NP);

for m=1:NT
    for i=1:NF
        for j=1:NW
            for k =1:NP
                
                Cpg = FCP_G.eval([PINL(m,i,j,k) Tin(m)]);
                Cpo = FCP_L.eval([PINL(m,i,j,k) Tin(m)]);
                
                WG = Fin(i)*Win(j);
                WO = (1-Fin(i))*Win(j);
                
                Hin_est(m,i,j,k) = T1(m,i,j,k)*(Cpg*WG + Cpo*WO);
                
                Hin_est(m,i,j,k) = Saturate(Hin_est(m,i,j,k),H_grid(1),H_grid(end));

                Error(m,i,j,k) = TTOP(m,i,j,k) - FTH_TOP.eval([Tin(m) Hin_est(m,i,j,k) Pout(k)]);
                %Error(m,i,j,k) = T2(m,i,j,k) - FT_TOP.eval([Tin(m) Fin(i) Win(j) Pout(k)]);
                
            end
            
        end
    end
        figure( 'Name', ['T2 Error for Tin=' num2str(Tin(m))] )
        clf
        surf(Win,Pout,reshape(Error(m,4,:,:),NW,NP)','FaceFin',0.3)
        xlabel('Win')
        ylabel('Pout')
        zlabel('T2 Error')
end
disp('Test completed.')
toc