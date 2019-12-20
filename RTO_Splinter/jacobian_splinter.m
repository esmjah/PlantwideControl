function JACx = jacobian_splinter( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L)
%COSTFUNC Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
x0 = x;
Cx0 = NonlinConScaled_splinter( x0,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L);

dx = zeros(n,1);
global espsilon

A1 = zeros(2,n);

for i = 1:n
    dx(i) = espsilon;
    Cx1 = NonlinConScaled_splinter( x0+dx,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L);
    A1(:,i) = (Cx1-Cx0)/dx(i);
    dx(i) = 0;
end
    
A2 = zeros(2,n);
for i = 1:n
    dx(i) = -espsilon;   
    Cx2 = NonlinConScaled_splinter( x0+dx,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L);
    A2(:,i) = (Cx2-Cx0)/dx(i);
    dx(i) = 0;
end

JACx = sparse((A1+A2)/2);

end

