% UTILDERIV Calculates derivative of utility wrt y.
% The root of this function represents ystar, the optimal choice of y. It
% also returns the second derivative, used by the NCPSOLVE function (from
% the CompEcon toolbox) to locate the root quickly. 
% 
% INPUTS
%   y       : earnings
%   lambda  : agent's laissez faire earnings
%   a       : lump sum tax subsidy
%   b       : marginal tax rate
% 
% OUTPUTS
%   uVal    : dU/dY at point y, equals zero at optimal y
%   uJac    : 2nd derivative (Jacobian of uVal) of utility at point y

function [uVal, uJac] = utilderiv(y,lambda,a,b)

global GAMMA SIGMA;

uVal = b*lambda.^(SIGMA + GAMMA - 1).*(a+b*y).^(-GAMMA) - y.^(SIGMA-1);
uJac = -GAMMA*b^2*lambda.^(SIGMA + GAMMA - 1).*(a+b*y).^(-GAMMA-1) - ...
       (SIGMA-1)*y.^(SIGMA-2);
