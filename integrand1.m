function [ y ] = integrand1( x )

y = exp(-(x-0.125).^2/0.125).*(1+cos(x*100)*0.1);

end

