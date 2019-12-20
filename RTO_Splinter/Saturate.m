function y = Saturate(x,xMin,xMax)
%SATURATE Summary of this function goes here
%   Detailed explanation goes here

y = x;

if(x>xMax)
    y = xMax;
end

if(x<xMin)
    y = xMin;
end


end

