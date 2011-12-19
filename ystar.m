% YSTAR Returns vector of agents' chosen earnings, given tax regime.
% Solution is computed using NCPSOLVE from the CompEcon toolbox, which
% finds the root of the derivative of u, subject to the constraint that y
% be nonnegative. (This is a complemetarity problem.)
% 
% INPUTS
%   lambdaArray : vector of agents' laissez faire earnings
%   a           : lump sum subsidy
%   b           : marginal tax rate
% 
% OUTPUTS
%   yStarArray  : vector of optimal earnings, given taxes
% 
% REQUIRED FUNCTIONS
%   UTILDERIV
%   NCPSOLVE

function yStarArray = ystar(lambdaArray,a,b)

global GAMMA SIGMA;

nAgents = size(lambdaArray,1);

y = lambdaArray(1);
y = lambdaArray(1).^(1/(SIGMA+GAMMA-1)); % guess laissez faire first
yStarArray = zeros(nAgents,1);

for i=1:nAgents
    lambda = lambdaArray(i);
    y = ncpsolve('utilderiv',0,inf,y,lambda,a,b); % requires y >= 0
    yStarArray(i) = y;
end 