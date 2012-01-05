% LAGRANGIANPHIBAR Returns vector derivatives of the Lagrangian. 
% The root of this function solves the optimal linear tax problem. Here we
% are interested heterogeneous phi, when phi is drawn from a normal
% distribution. We want to answer the following question: Given a vector of
% var(phi) (with mean zero) what mean(phi) is required to generate a given
% optimal marginal tax rate.
% 
% INPUTS
%   sol                 : vector [meanPhi; a; q], where a is the lump sum
%                         subsidy and q is the shadow social value of lump 
%                         sum subsidy. 
%   mtrTarget           : the target marginal tax rate
%   demeanedPhiArray    : variation from mean of phi_i distribution, so 
%                         that phiArray = meanPhi+varPhiArray.
%   lambdaArray         : vector of laissez faire earnings
% 
% OUTPUTS
%   fval                : Jacobian of derivative of Lagrangian wrt sol -- 
%                         this is [0; 0; 0] at solution. 
% 
% REQUIRED FUNCTIONS
%   YSTAR
%   YJACOBIAN

function fval = lagrangianphibar(sol,mtrTarget,demeanedPhiArray,lambdaArray)

global GAMMA SIGMA;

meanPhi = sol(1);
a = sol(2);
q = sol(3);

b = 1-mtrTarget;

phiArray = meanPhi+demeanedPhiArray;
thetaArray = lambdaArray.^((SIGMA+GAMMA-1)*phiArray);

yStarArray = ystar(lambdaArray,a,b);    % agents' earning choices given tax

mu = (SIGMA-1)/(SIGMA+GAMMA-1);         % for pure preferences

dUdA = thetaArray.^(1-mu).*(a+b*yStarArray).^(-GAMMA);
dUdB = yStarArray.*dUdA;

yJac = yjacobian(lambdaArray,yStarArray,a,b);
dYdA = yJac(:,1);
dYdB = yJac(:,2);

dLdA = sum(dUdA - q*(1-(1-b)*dYdA));
dLdB = sum(dUdB - q*(yStarArray - (1-b)*(dYdB)));

budgetConstraint = -a + (1-b)*mean(yStarArray);

fval = [dLdA; dLdB; budgetConstraint];
