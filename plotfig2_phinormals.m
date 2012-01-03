% OPTIMAL LINEAR TAX WITH HETEROGENEOUS PHI DRAWN FROM MEAN-1 NORMAL DIST 
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
% theta_i = lambda_i^((sigma+gamma-1)*(1-phi_i)) and 
% w_i = lambda_i^((sigma+gamma-1)*phi_i/sigma).
% 
% This script finds the optimal linear tax regime (described by c=a+b*y)
% when phi is drawn from a normal distribution with mean 1 and variance in
% [0,0.5]. By "optimal", we mean the tax regime which maximizes utilitarian
% social welfare over the pure preferences specification above. (Note that
% other, non-pure-preference specifications can also be found by changing
% the value of mu in LAGRANGIAN.)
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
loadcompecon();

% Customizeable options:
nAgents = 500;

global GAMMA SIGMA;     % declare global parameters
GAMMA = 1;
SIGMA = 3;

% Generate matrix of phi distributions
meanPhi = 1;
varPhiArray = 0:0.05:0.5;
nPhiDists = size(varPhiArray,2);
meanPhiMatrix = meanPhi.*ones(nAgents,nPhiDists);
varPhiMatrix = normrnd(0,1,nAgents,1)*varPhiArray;
phiMatrix = meanPhi + varPhiMatrix;

lambdaArray = simulateagents(nAgents);      % simulate agents
options = optimset('Display', 'iter');       % set optimization options

sol = [0; 1; 0.01]; % first guess of solution to lagrangian: [a; b; q]
mtrArray = zeros(nPhiDists,1);
for iPhi=1:nPhiDists
    phiArray = phiMatrix(:,iPhi);
    disp(['iteration ' num2str(iPhi)]);             % display iter no.
    thetaArray = lambdaArray.^((SIGMA + GAMMA - 1)*(1-phiArray));
    
    % Find root of Lagrangian derivative
    sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);
    mtrArray(iPhi) = 1-sol(2);                  % mtr = 1-b
end

% Plot results
line(varPhiArray,mtrArray,'LineWidth',2,'Color',[0 0 1]);%figure(gcf)
axis([0 0.5 0 1]);
xlabel('var(\phi)');
ylabel('Marginal tax rate');
