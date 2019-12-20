function [FP_INL,FT_INL,FT_TOP] = Build_pipeline_functions()
%BUILD_PIPELINE_FUNCTIONS Summary of this function goes here
%   Detailed explanation goes here

%load PipelineData.mat
load Pipeline_Olga2014_Data.mat


NT = length(Tin);
NF = length(Fin);
NW = length(Win);
NP = length(Pout);


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
                                
                PINLs(n) = PINL(m,i,j,k);
                TINLs(n) = TINL(m,i,j,k);
                TTOPs(n) = TTOP(m,i,j,k);
                
                n = n + 1;
            end
            
        end
    end
end

FP_INL = BSplineBuilder(xs,PINLs, 3).build();
disp('FP is done')
FT_INL = BSplineBuilder(xs,TINLs, 3).build();
disp('FTINL is done')
FT_TOP = BSplineBuilder(xs,TTOPs, 3).build();
disp('FTOP is done')

end

