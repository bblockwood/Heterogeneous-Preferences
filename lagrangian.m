% LAGRANGIAN Returns vector derivatives of the Lagrangian.
% The root of this function solves the optimal linear tax problem. 
% 
% INPUTS
%   sol         : vector [a; b; q], where post-tax income is a+b*y, and q
%                 is the shadow social value of lump sum subsidy. 
%   thetaArray  : vector of taste parameters
%   lambdaArray : vector of lambda parameters (laissez faire earnings)
% 
% OUTPUTS
%   fval        : Jacobian of derivative of Lagrangian wrt sol -- this
%                 is [0; 0; 0] at solution. 
% 
% REQUIRED FUNCTIONS
%   YSTAR
%   YJACOBIAN

function fval = lagrangian(sol,thetaArray,lambdaArray)

global GAMMA SIGMA;

a = sol(1);
b = sol(2);
q = sol(3);

yStarArray = ystar(lambdaArray,a,b);    % agents' earning choices given tax

mu = (SIGMA-1)/(SIGMA+GAMMA-1);         % for pure preferences

% dUdA = thetaArray.^(1-mu).*(a+b*yStarArray).^(-GAMMA);
dUdA = thetaArray.^(1-mu).*(a+b*yStarArray).^(-GAMMA);

dUdB = yStarArray.*dUdA;

yJac = yjacobian(lambdaArray,yStarArray,a,b);
dYdA = yJac(:,1);
dYdB = yJac(:,2);

dLdA = sum(dUdA - q*(1-(1-b)*dYdA));
dLdB = sum(dUdB - q*(yStarArray - (1-b)*(dYdB)));

budgetConstraint = -a + (1-b)*mean(yStarArray);

fval = [dLdA; dLdB; budgetConstraint];
