N = 100;
sigmas =(0.05:0.05:2.5)*pi/4;
M = length(sigmas);
var0 = zeros(1, M);
var1 = zeros(1, M);
var2 = zeros(1, M);
r0 = zeros(1, N);
r1 = zeros(1, N);
r2 = zeros(1, N);

idx = 1;
for sigma = sigmas,
for k=1:N,
    
% Parameters for the model setup

% The model scene:
% We define a light source at y=1, ranging from x=0 to x=1. The light
% source emits radiance of 1. We define a BRDF as a Gaussian, centered
% around direction at angle pi/4. This models a glossy surface seen from a
% viewing direction at angle -pi/4. The standard deviation of the Gaussian
% models surface roughness. The integral we compute is simply 
% int[-pi/2,pi/2] BRDF(w) L(w) dw. 
%
% Note: the light geometry implies L(w)=1 for 0<w<pi/4, and 0 otherwise.

% Number of samples
n = 500;

% Parameters of BRDF
mu = pi/4; % Mean, you should not change this

% Sigma corresponds to surface roughness in our experiment here. Use
% different sigma values to experiment with varying surface roughness.
%sigma = pi/4 * 0.25;


%%
% Sample BRDF

% Sample over directions with Gauss distribution
% Samples (directions)
w = gaussSamples( mu , sigma, n, -pi/2, pi/2, false);
% Their probability densities
p = gaussSampleProbabilities( w, mu, sigma , -pi/2, pi/2);
prob0 = @(x)gaussSampleProbabilities( x, mu, sigma , -pi/2, pi/2);
p0 = p;
w0 =w;
% BRDF values
brdf = p;

% The light is defined to be in the range of angles [0,pi/4]
% Since the pdf of the samples corresponds to the sample values, after
% dividing sample values by pdfs we obtain 1s!
f = (w > 0 & w < pi/4);
f0 = f;
% Monte Carlo estimate
r = sum(f)/n;
r0(k) = r;
%%
% Sample area light

% Random samples on light source. Light source is defined at at y=1, going
% from x=0 to x=1. We choose random x locations in the range [0,1] for the
% light samples
s = rand(1,n);

% Corresponding directions (angles)
w = atan(s);
w1=w;
% Densities over directions
% This implements the conversion from Equation 8.10 in the thesis. Because 
% we are in 2D (instead of 3D), there is a factor r instead of r^2.
% We compute p(w) = r/cos * p(x), note that p(x)=1 in our case.
p = sqrt(1+s.^2)./(1./sqrt(1+s.^2));
prob1 = @(x)sqrt(1+tan(x).^2)./(1./sqrt(1+tan(x).^2));

p1 = p;
% BRDF values
brdf = gaussSampleProbabilities(w, mu, sigma, -pi/2, pi/2);

% Sample values, BRDF divided by probability densities
f = brdf./p;
f1 = f;
% Monte Carlo estimate
r = sum(f)/n;
r1(k) = r;
%%

w0_hat = (p0*n) ./(n*prob0(w0)+n*prob1(w0)) ;
w1_hat = (p1*n) ./(n*prob0(w1)+n*prob1(w1));

r = (1/n)*sum(f0.*w0_hat + f1.*w1_hat);
r2(k) = r;
%%
% Visualization
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
hold on
line([-1 0],[1 0], 'Color','b')
line([0 1],[1 1], 'Color','c')
plot(sin(w),cos(w),'.');
plot(sin(w).*f,cos(w).*f,'.red')
plot(sin(w).*brdf,cos(w).*brdf,'.green')
axis([0 2 0 2]);
axis equal
legend('viewing direction','light source','sample directions','sample values (shown as distance to origin)','brdf values (shown as distance to origin)')
t = sprintf('Monte Carlo estimate r=%f', r);
title(t)
%%

% just plotting weights for last sigma in above's most outer iteration
figure
title('plot against p0')
plot(p0,w0_hat,'.',p0,w1_hat, '.')
legend('w0_hat', 'w1_hat')

figure
plot(p1,w0_hat,'.',p1,w1_hat, '.')
legend('w0_hat', 'w1_hat')