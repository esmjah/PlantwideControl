function JACx = jacobian( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FH_G,FH_L,FTH_TOP)
%COSTFUNC Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
x0 = x;
Cx0 = NonlinConScaled( x0,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FH_G,FH_L,FTH_TOP);

dx = zeros(n,1);
eps = 1e-6;

A1 = zeros(2,n);

for i = 1:n
    dx(i) = eps;
    Cx1 = NonlinConScaled( x0+dx,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FH_G,FH_L,FTH_TOP);
    A1(:,i) = (Cx1-Cx0)/eps;
    dx(i) = 0;
end
    
eps = -eps;
A2 = zeros(2,n);
for i = 1:n
    dx(i) = eps;   
    Cx2 = NonlinConScaled( x0+dx,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FH_G,FH_L,FTH_TOP);
    A2(:,i) = (Cx2-Cx0)/eps;
    dx(i) = 0;
end

JACx = sparse((A1+A2)/2); %sparse(A1);

end

