function [ s ] = gaussSamples( mu , sigma, n, i_min, i_max, stratify)
%GAUSSSAMPLES mu sigma describe the normal distribution, n the number of
%samples, i_min, i_max, the target interval


    %gauss = @(x) 1/(sigma*sqrt(2*pi))*exp(-(x-mu).^2/(2*sigma^2)) ;
    s_to_u = @(s) (erf((s-mu)/(sigma*sqrt(2)))+1)/2;
    %u_to_s = @(u) mu + sigma*sqrt(2)*erfinv(2*u-1);
    
    if(stratify)
        u = ((0:n-1)+ rand(1,n))./n;
    else
        u = rand(1,n);
    end
    
    u_min = s_to_u(i_min);
    u_max = s_to_u(i_max);
    u = u_min + u *(u_max-u_min);


    % Note: the quantile f, that is the probability that a random sample of
    % a distribution is smaller than f, of the normal distribution is
    % mu + sigma*sqrt(2)*erfinv(2*f-1)

    % Samples distributed according to 1/sqrt(pi)*exp(-x.^2), that is, 
    % a normal distribution with mu=0 and sigma=1/sqrt(2)
    s = mu + sigma*sqrt(2)*erfinv(2*u-1); % inverse of cfd / Quatntile

end

