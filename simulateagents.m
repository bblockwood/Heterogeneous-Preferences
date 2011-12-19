% SIMULATEAGENTS Simulates agents
% 
% INPUTS
%   nAgents     : whaddayathink?
% 
% OUTPUTS
%   ylfArray    : vector of laissez-faire earnings

function ylfArray = simulateagents(nAgents)

mu = 1.6;
sd = 0.8;

% Simulate lambdas (calibrated to US inc distribution in 10,000's)
ylfArray = sort(lognrnd(mu, sd, nAgents, 1));