% FINDYJACROW Computes one row of the the Jacobian of yStarArray. 
% This is a vector of equations implicitly characterizing the Jacobian of a
% single agent's yStar wrt the vector describing the tax regime, [a; b].
% The root of this function is the actual Jacobian, so that when result =
% [0; 0] then x = [dYstar/da; dYstar/db].
% 
% INPUTS
%   x       : test value for [dYstar/da; dYstar/db]
%   yStar   : agent's chosen earnings given tax regime [a b]
%   lambda  : agent's laissez faire earnings
%   a       : lump sum tax subsidy
%   b       : marginal tax rate
% 
% OUTPUTS
%   result  : implicit equations characterizing Jacobian. Equals [0; 0]
%             when x is the correct [dYstar/da; dYstar/db].

function result = findyjacrow(x,yStar,lambda,a,b)

global GAMMA SIGMA;

dyda = x(1);
dydb = x(2);

rootA = -GAMMA*b*lambda.^(SIGMA + GAMMA - 1).*...
        (a+b*yStar).^(-GAMMA-1).*(1+b*dyda) - ...
        (SIGMA-1)*yStar.^(SIGMA-2).*dyda;

rootB = lambda.^(SIGMA + GAMMA - 1).*(a+b*yStar).^(-GAMMA) - ...
        GAMMA*b*lambda.^(SIGMA + GAMMA - 1).*...
        (a+b*yStar).^(-GAMMA-1).*(yStar+b*dydb) - ...
        (SIGMA-1)*yStar.^(SIGMA-2).*dydb;

result = [rootA; rootB];
