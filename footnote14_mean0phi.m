% OPTIMAL LINEAR TAX WITH HETEROGENEOUS PHI DRAWN FROM MEAN-0 DISTRIBUTIONS
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
% This script finds the optimal linear tax regime (described by c=a+b*y)
% when phi is drawn from different mean-0 distributions. The goal is to
% find out whether heterogeneous phi tends to yield a more redistributive
% optimal tax regime than constant phi, regardless of the shape of the
% distribution. By "optimal", we mean the tax regime which maximizes
% utilitarian social welfare over the preference neutral cardinalization.
% (Note that other, non-pure-preference specifications can also be found by
% changing the value of mu in LAGRANGIAN.)
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
%   SIMULATEAGENTS
%   UTILDERIV
%   YJACOBIAN
%   YSTAR

clear all;

% Customizeable options:
nAgents = 1000;

global GAMMA SIGMA;     % declare global parameters
GAMMA = 1;
SIGMA = 3;

lambdaArray = simulateagents(nAgents);      % simulate agents
options = optimset('Display', 'iter');       % set optimization options

% Simulate constant phi=0
phi=0;
thetaArray = lambdaArray.^((SIGMA + GAMMA - 1)*phi); % (vector of ones)
sol = [2.6; 0.5; 0.18]; % first guess of solution to lagrangian: [a; b; q]
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);
disp(['optimal mtr for constant phi=0: ' num2str(1-sol(2))]);


% Lognormal distributions
meanNormal = 0; % params for normal dist associated w/ logN dist
sdNormal = 0.5;
meanLogN = exp(meanNormal+sdNormal^2/2);

% Simulate, X ~ logN(0,1), with phi=X-mean(X) (right skew)
% (Helps to run this section with smaller # of agents, ~100, to see t<0.)
phiArray = lognrnd(meanNormal,sdNormal,nAgents,1)-meanLogN;
thetaArray = lambdaArray.^((SIGMA + GAMMA - 1)*phiArray);
sol = [-1.7; 1.6; 0.6]; % first guess of solution to lagrangian: [a; b; q]
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options)
disp(['optimal mtr for phi=X-E(X), X~LogN(0,0.5): ' num2str(1-sol(2))]);

% Simulate, X ~ logN(0,1), with phi=mean(X)-X (left skew)
phiArray = -lognrnd(meanNormal,sdNormal,nAgents,1)+meanLogN;
thetaArray = lambdaArray.^((SIGMA + GAMMA - 1)*phiArray);
sol = [2; 0.6; 0.2]; % first guess of solution to lagrangian: [a; b; q]
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options)
disp(['optimal mtr for phi=E(X)-X, X~LogN(0,0.5): ' num2str(1-sol(2))]);
