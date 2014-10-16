%% 
% Generate samples distributed according to a normal distribution with
% sigma = 1/sqrt(2). The code here uses the inversion method (Sec. 2.4.2 in
% the thesis) to obtain the desired distribution.

% Number of samples
n = 100000;

% Uniformly distributed samples in [-1,1]
u = 2*rand(1,n)-1;

% Note: the quantile f, that is the probability that a random sample of
% a distribution is smaller than f, of the normal distribution is
% mu + sigma*sqrt(2)*erfinv(2*f-1)

% Samples distributed according to 1/sqrt(pi)*exp(-x.^2), that is, 
% a normal distribution with mu=0 and sigma=1/sqrt(2)
s = erfinv(u);

% Visualize
[h x] = hist(s, n/1000);
dx = diff(x(1:2));
bar(x,h/sum(h*dx));
hold on
plot(x, 1/sqrt(pi)*exp(-x.^2),'red')