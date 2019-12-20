function Gx = gradient_splinter(x,Lambda,WL_w1,WL_w2)
%COSTFUNC Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
x0 = x;
dx = zeros(n,1);
global espsilon

G1 = zeros(n,1);
Jx0 = CostFuncScaled_splinter(x0,Lambda,WL_w1,WL_w2);


for i = 1:n
    dx(i) = espsilon;
    Jx1 = CostFuncScaled_splinter(x0+dx,Lambda,WL_w1,WL_w2);
    G1(i) = (Jx1-Jx0)/dx(i);
    dx(i) = 0;
end
    
G2 = zeros(n,1);
for i = 1:n
    dx(i) = -espsilon;
    Jx1 = CostFuncScaled_splinter(x0+dx,Lambda,WL_w1,WL_w2);
    G2(i) = (Jx1-Jx0)/dx(i);
    dx(i) = 0;
end

Gx = G1; %(G1+G2)/2 %G1;

end

