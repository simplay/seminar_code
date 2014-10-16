function [ res, weights, probs ] = integrate( sigma, n )
res = zeros(1,3);
weights = zeros(1,n,2);
probs = zeros(1,n,2);

%INTEGRATE Summary of this function goes here
%   Detailed explanation goes here

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
probs(:,:,1) = p0;
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
res(1) = r;
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
probs(:,:,2) = p1;
% BRDF values
brdf = gaussSampleProbabilities(w, mu, sigma, -pi/2, pi/2);

% Sample values, BRDF divided by probability densities
f = brdf./p;
f1 = f;
% Monte Carlo estimate
r = sum(f)/n;
res(2) = r;
%%

w0_hat = (p0*n) ./(n*prob0(w0)+n*prob1(w0)) ;
w1_hat = (p1*n) ./(n*prob0(w1)+n*prob1(w1));

r = (1/n)*sum(f0.*w0_hat + f1.*w1_hat);
res(3) = r;

weights(:,:,1) = w0_hat;
weights(:,:,2) = w1_hat;

end

