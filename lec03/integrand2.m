function [ y ] = integrand2( x )

y = exp(-(x-0.125).^2/0.125).*(1+cos(x*100)*0.1)+ 3*exp(-(x-0.8).^2/0.01).*(1+cos(x*200)*0.2);

end

