function [ y ] = sample( n )
%SAMPLE Summary of this function goes here
%   Detailed explanation goes here

x = rand(n, 1);
y = integrand1(x);
y = var(x,y);
end

