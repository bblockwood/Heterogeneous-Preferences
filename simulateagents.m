% SIMULATEAGENTS Simulates agents
% 
% INPUTS
%   nAgents     : whaddayathink?
% 
% OUTPUTS
%   ylfArray    : vector of laissez-faire earnings

function ylfArray = simulateagents(nAgents)

mu = 1.65;
sd = 0.75;

% Simulate lambdas (calibrated to US inc distribution in 10,000's)
rng(1);
ylfArray = sort(lognrnd(mu, sd, nAgents, 1));
