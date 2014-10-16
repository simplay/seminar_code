function [ p ] = gaussSampleProbabilities( x, mu, sigma , i_min, i_max)
%GAUSSSAMPLEPROBABILITIES Summary of this function goes here
%   Detailed explanation goes here
gauss = @(x) 1/(sigma*sqrt(2*pi))*exp(-(x-mu).^2/(2*sigma^2));
p = gauss(x)/(normcdf(i_max,mu,sigma)-normcdf(i_min,mu,sigma));

end

