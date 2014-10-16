function [ t2] = sample2( n )
%SAMPLE2 Summary of this function goes here
%   Detailed explanation goes here
    x = intervalSample(n);
    y = integrand1(x);
    t1 = sum(y)/n;
    
    t2 = var(x,y);

end

