function [ output_args ] = var( x, y )
%X Summary of this function goes here
%   Detailed explanation goes here
    output_args = (x.^2)'*y-(x'*y)^2;
    
end

