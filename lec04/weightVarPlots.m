function weightVarPlots( n, N, sigmaCount )
%WEIGHTVARPLOTS Plot How error decreases for different strats + weights
%   n number of samples used for integration
%   N number of samples for computing the variance
%   sigmaCount number of partitions splitting the sigma range [0.05, 2.5]
%
% E.g. weightVarPlots(500, 400, 10)

sigmas = linspace(0.05,2.5,sigmaCount)*pi/4;
M = length(sigmas);
var0 = zeros(1, M);
var1 = zeros(1, M);
var2 = zeros(1, M);
r0 = zeros(1, N);
r1 = zeros(1, N);
r2 = zeros(1, N);

w0_hat = zeros(1,n,2); w1_hat = zeros(1,n,2); p0 = zeros(1,n,2); p1 = zeros(1,n,2);

idx = 1;
for sigma = sigmas,
    for k=1:N,
        [res, weights, probs] = integrate(sigma, n, 0);
        r0(k) = res(1); r1(k) = res(2); r2(k) = res(3);
        w0_hat = weights(:,:,1); w1_hat = weights(:,:,2);
        p0 = probs(:,:,1); p1 = probs(:,:,2);
    end
    var0(idx) = var(r0);
    var1(idx) = var(r1);
    var2(idx) = var(r2);
    
    idx = idx + 1;
end

figure
hold on
plot(sigmas,var0,'b', sigmas,var1,'k', sigmas,var2,'r')
legend('brd', 'light', 'comb')
hold off

figure
title('plot against p0')
plot(p0,w0_hat,'.',p0,w1_hat, '.')
legend('w0_hat', 'w1_hat')

figure
plot(p1,w0_hat,'.',p1,w1_hat, '.')
legend('w0_hat', 'w1_hat')
    
end

