% PLOTS OF ISOQUANTS OF MTR AND OF VAR(LN(MRS))
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
% corresponding phibar which would generate the same tax regime. We plot
% isoquants in mean(phi), var(phi) space for both the marginal tax rate and
% variance of reported MRS (theta). 
% 
% For more on the derivation of the formulas herein, see the accompanying
% file notes_numeric_optimization.pdf.
% 
% REQUIRED PACKAGES
%   COMPECON (www4.ncsu.edu/~pfackler/compecon/toolbox.html)
% 
% REQUIRED FUNCTIONS
%   FINDYJACROW
%   LAGRANGIAN
%   LOADCOMPECON
%   SIMULATEAGENTS
%   UTILDERIV
%   YJACOBIAN
%   YSTAR

clear all;

% Customizeable options:
nAgents = 500;
global GAMMA SIGMA;
GAMMA = 1;
SIGMA = 3;

phibar = 0.2;                           % constant phi of interest
options = optimset('Display', 'iter','MaxFunEvals',1000);

lambdaArray = simulateagents(nAgents);  % simulate agents
thetaArray = lambdaArray.^((SIGMA + GAMMA - 1)*phibar);

% Generate matrix of phi deviations
varPhiArray = [0.01 0.05:0.05:0.6];
nPhiDists = size(varPhiArray,2);
variationPhi = normrnd(0,1,nAgents,1);

% Find MTR isoquants
sol = [0.5; 1; 0.5]; % first guess of solution to lagrangian: [a; b; q]
sol = fsolve(@(x)lagrangian(x,thetaArray,lambdaArray),sol,options);

mtrTarget = 1-sol(2);   % optimal mtr given constant phibar
aGuess = sol(1);
qGuess = sol(3);

mtrIso = zeros(nPhiDists,1);        % we're solving for this
sol = [phibar; aGuess; qGuess];     % initial guess is phibar solution
varLnThetaArray = zeros(nPhiDists,1);
corrArray = zeros(nPhiDists,1);

for iPhiDist = 1:nPhiDists
    varPhi = varPhiArray(iPhiDist);
    disp(['var(phi) = ' num2str(varPhi)]);  % echo status   
    demeanedPhiArray = variationPhi*varPhi;
       
    if varPhi~=0    % if varPhi = 0, mtr is just the precalculated case
        sol = fsolve(@(x)...
            lagrangianphibar(x,mtrTarget,demeanedPhiArray,lambdaArray),...
            sol,options);
    end 

    mtrIso(iPhiDist) = sol(1);                      % store mean phi
    earnings = ystar(lambdaArray,sol(2),sol(3));    % earnings given taxes
    
    % For this MTR, find var(ln(theta)) and corr(theta,lambda)
    phiArray = sol(1) + demeanedPhiArray;
    thetaArray = lambdaArray.^((SIGMA+GAMMA-1)*phiArray);
    corrMrsEarnings = corr([thetaArray earnings]);
    varLnThetaArray(iPhiDist) = var(log(thetaArray));
    corrArray(iPhiDist) = corrMrsEarnings(2,1);    
end

% Plot constant mtr in var(ln(mrs)) & corr(mrs,earnings) space
line(varLnThetaArray,corrArray,'LineWidth',2,'Color',[0 0 1]);
axis([0 1 0 1]);
xlabel('var(ln(\theta))');
ylabel('corr(\theta,y)');

% Add text boxes
str1 = {'Optimal \tau < \tau*'};
text(0.6,0.75,str1,'FontSize',14);
str2 = {'Optimal \tau > \tau*'};
text(0.3,0.4,str2,'FontSize',14);

% display tax rate
disp(mtrTarget);
