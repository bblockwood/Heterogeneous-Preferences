% PLOTS HISTOGRAMS OF WAGE DISTRIBUTIONS FOR VARIOUS PHI DISTRIBUTIONS
% 
% This script shows how preference heterogeneity affects the optimal
% marginal tax rate (mtr). We suppose agents have the utility function
% U_i(c,l) = theta_i*u(c)-v(l), where u(c) = (c^(1-gamma)-1)/(1-gamma),
% v(l) = l^sigma/sigma, and l = y/w_i; thus there are two dimensions of
% heterogeneity: skill (w_i) and tastes (theta_i). This utility function is
% observationally equivalent to the pure preference cardinalization
% U_i(c,l) = theta_i^(1-mu)*u(c)-theta^(-mu)*v(l), with
% mu=(sigma-1)/(sigma+gamma-1). Agents are characterized by
% lambda_i=(theta_i*w_i^sigma)^(1/(sigma+gamma-1)), which is their optimal
% income given a laissez faire tax regime, with 
% theta_i = lambda_i^((sigma+gamma-1)*phi_i) and 
% w_i = lambda_i^((sigma+gamma-1)*(1-phi_i)/sigma).
% 
% This script finds plots the frequencies of w for various phi
% distributions, and calculates the optimal marginal tax rate associated
% with each distribution. 
% 
% For more on the derivation of the formulas herein, see the accompanying
% file notes_numeric_optimization.pdf.
% 
% REQUIRED PACKAGES
%   COMPECON (www4.ncsu.edu/~pfackler/compecon/toolbox.html)
% 
% REQUIRED FUNCTIONS
%   COMPUTEYJACROW
%   LAGRANGIAN
%   LOADCOMPECON
%   SIMULATEAGENTS
%   UTILDERIV
%   YJACOBIAN
%   YSTAR

clear all;
clc;

% Customizeable options:
nAgents = 1000;
nPhi = 50;
nAgents = 50;       % for optimal tax calculation
nPhi = 10;          % for optimal tax calculation
nBins = 40;

global GAMMA SIGMA;     % declare global parameters
GAMMA = 1;
SIGMA = 3;

% Draw lambdas for constant phi
mu = 1.65;
sd = 0.75;
stepsize = 1/(nAgents*nPhi);
draws = (stepsize/2:stepsize:1-stepsize/2)';
lambdaArray = logninv(draws,mu,sd);

wUB = 15;                                   % upper bound of horiz axis
% freqUB = 10000;
axisBnds = [0 wUB 0 freqUB];
options = optimset('Display', 'iter');      % set optimization options


%% Mirrlees benchmark: E(phi) = 0, var(phi) = 0
phiBar = 0;
wArray = lambdaArray.^((SIGMA+GAMMA-1)*(1-phiBar)/SIGMA);

subplot(2,2,1);
n = hist(wArray(wArray < wUB),nBins);
freqUB = 1.5*max(n);
axis(axisBnds);
title('E(\phi) = 0, var(\phi) = 0');

% Find optimal MTR
sol = [2.8; 0.5; 0.2]; % first guess of solution to lagrangian: [a; b; q]
thetaArray = lambdaArray.^((SIGMA + GAMMA - 1)*phiBar);
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);
mtrArray(1) = 1-sol(2)


%% E(phi) = 0.5, var(phi) = 0
phiBar = 0.5;
wArray = lambdaArray.^((SIGMA+GAMMA-1)*(1-phiBar)/SIGMA);

subplot(2,2,2);
hist(wArray(wArray < wUB),nBins);
axis(axisBnds);
title('E(\phi) = 0.5, var(\phi) = 0');

% Find optimal MTR
sol = [1.8; 0.7; 0.4]; % first guess of solution to lagrangian: [a; b; q]
thetaArray = lambdaArray.^((SIGMA + GAMMA - 1)*phiBar);
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);
mtrArray(2) = 1-sol(2)


%% Prep for variable phi
stepsize = 1/nAgents;
draws = (stepsize/2:stepsize:1-stepsize/2)';
lambdaArray = logninv(draws,mu,sd);

% Generate normal draws for phi distribution
stepsize = 1/nPhi;
draws = (stepsize/2:stepsize:1-stepsize/2)';
normArray = norminv(draws,0,1);


%% E(phi) = 0, var(phi) = 0.25
lambdaMat = repmat(lambdaArray,1,nPhi);
phiArray = 0.5*normArray;
phiMat = repmat(phiArray',nAgents,1);
wMat = lambdaMat.^((SIGMA+GAMMA-1)*(1-phiMat)/SIGMA);
wArray = wMat(:);
subplot(2,2,3);
hist(wArray(wArray < wUB),nBins);
axis(axisBnds);
title('E(\phi) = 0, var(\phi) = 0.25');

% Find optimal MTR
sol = [2.2; 0.6; 0.2]; % first guess of solution to lagrangian: [a; b; q]
thetaMat = lambdaMat.^((SIGMA + GAMMA - 1)*phiMat);
thetaArray = thetaMat(:);
lambdaArray = lambdaMat(:);
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);
mtrArray(3) = 1-sol(2)



%% E(phi) = 0.5, var(phi) = 0.25
lambdaMat = repmat(lambdaArray,1,nPhi);
phiArray = 0.5*normArray+0.5;
phiMat = repmat(phiArray',nAgents,1);
wMat = lambdaMat.^((SIGMA+GAMMA-1)*(1-phiMat)/SIGMA);
wArray = wMat(:);
subplot(2,2,4);
hist(wArray(wArray < wUB),nBins);
axis(axisBnds);
title('E(\phi) = 0.5, var(\phi) = 0.25');

% Find optimal MTR
sol = [1.8; 0.7; 0.4]; % first guess of solution to lagrangian: [a; b; q]
thetaMat = lambdaMat.^((SIGMA + GAMMA - 1)*phiMat);
thetaArray = thetaMat(:);
lambdaArray = lambdaMat(:);
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);
mtrArray(4) = 1-sol(2)
