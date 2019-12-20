function Gx = gradient(x,Lambda,WL_w1,WL_w2)
%COSTFUNC Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
x0 = x;
dx = zeros(n,1);
eps = 1e-6;

G1 = zeros(n,1);
Jx0 = CostFuncScaled(x0,Lambda,WL_w1,WL_w2);


for i = 1:n
    dx(i) = eps;
    Jx1 = CostFuncScaled(x0+dx,Lambda,WL_w1,WL_w2);
    G1(i) = (Jx1-Jx0)/eps;
    dx(i) = 0;
end
    
eps = -eps;
G2 = zeros(n,1);
for i = 1:n
    dx(i) = eps;
    Jx1 = CostFuncScaled(x0+dx,Lambda,WL_w1,WL_w2);
    G2(i) = (Jx1-Jx0)/eps;
    dx(i) = 0;
end

Gx = (G1+G2)/2; %G1;

end

