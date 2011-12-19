% YJACOBIAN Computes the Jacobian of the yStar vector wrt tax regime [a b].
% Solution is computed using BROYDEN from the CompEcon toolbox. 
% 
% INPUTS
%   lambdaArray : vector of agents' laissez faire earnings
%   yStarArray  : vector of agents' optimal earnings
%   a           : lump sum subsidy
%   b           : marginal tax rate
% 
% OUTPUTS
%   yJac        : Jacobian matrix, dim = (nAgents, 2)
% 
% REQUIRED FUNCTIONS
%   FINDJACROW
%   BROYDEN

function yJac = yjacobian(lambdaArray,yStarArray,a,b)

nAgents = size(lambdaArray,1);

x = [0; 0]; % starting guess
yJac = zeros(nAgents,2);
for i=1:nAgents
    lambda = lambdaArray(i);
    yStar = yStarArray(i);
    x = broyden('findyjacrow',x,yStar,lambda,a,b); % x=[dy_i/da dy_i/db]
    yJac(i,:) = x';
end