% DISTRIBUTIONS OF PHI GENERATING SAME OPTIMAL MTR AS CONSTANT PHI=0
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
% In this script we allow phi_i to vary. We suppose phi_i is drawn from a
% normal distribution. For any such distribution of phi, there is a
% corresponding phibar which would generate the same tax regime. Here we
% ask what combinations of mean(phi) and var(phi) would generate the same
% optimal mtr as a given constant phibar. We generate a plot of the
% mean(phi) necessary to generate that mtr for a range of var(phi) between
% 0 and 0.5.
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

% Customizeable options:
nAgents = 500;
global GAMMA SIGMA;                     % declare global parameters
GAMMA = 1;
SIGMA = 3;
phibar = 0;                             % constant phi of interest
options = optimset('Display', 'off');   % set optimization options

lambdaArray = simulateagents(nAgents);                  % simulate agents
thetaArray = lambdaArray.^((SIGMA+GAMMA-1)*phibar);  % "tastes" vector

% Find optimal mtr given constant phibar
sol = [0; 1; 0.01]; % first guess of solution to lagrangian: [a; b; q]
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);
mtrTarget = 1-sol(2);
aGuess = sol(1);
qGuess = sol(3);

% Generate matrix of phi distributions
sdPhiArray = 0:0.1:0.5;
nPhiDists = size(sdPhiArray,2);
rng(1);
variationPhi = normrnd(0,1,nAgents,1);

meanPhiArray = zeros(nPhiDists,1);
sol = [phibar; aGuess; qGuess];     % initial guess is phibar solution
for iPhiDist = 1:nPhiDists
    sdPhi = sdPhiArray(iPhiDist);
    disp(['sd(phi) = ' num2str(sdPhi)]);   % echo status
    
    if sdPhi==0, meanPhiArray(iPhiDist)=phibar; continue, end % by design
    
    demeanedPhiArray = variationPhi*sdPhi;
    sol = fsolve(...
        @(x) lagrangianphibar(x,mtrTarget,demeanedPhiArray,lambdaArray),...
        sol,options);
    meanPhiArray(iPhiDist) = sol(1);
end

% Plot results
varPhiArray = sdPhiArray.^2;
line(varPhiArray,meanPhiArray,'LineWidth',2,'Color',[0 0 1]);
axis([0 0.25 -0.5 1]);
xlabel('var(\phi)');
ylabel('Required mean(\phi)');

% Add text boxes
str1(1) = {'Less redistribution'};
str1(2) = {'than conventional model'};
text(0.07,0.1,str1,'FontSize',14);
str2(1) = {'More redistribution'};
str2(2) = {'than conventional model'};
text(0.02,-0.3,str2,'FontSize',14);
