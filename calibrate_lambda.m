% CALIBRATE_LAMBDA Calibrates parameters of lambda distribution. 
% Finds the parameters of the lognormal distribution that minimize the
% distance (sum of squares) between US income quintiles and simulated data.

global GAMMA SIGMA;
GAMMA = 1;
SIGMA = 3;
nAgents = 5000;

% Lower bounds for deciles 2:5
% source: http://pubdb3.census.gov/macro/032005/hhinc/new05_000.htm
LBs = [1.85; 3.47; 5.53; 8.80];
prc = [20; 40; 60; 80];    % percentiles to be matched

% Starting values
mu0 = 1;
sd0 = 1;
th0 = [mu0; sd0];       % combined parameter vector

% Tax parameters, from Kotlikoff and Rapson, 2006, NBER WP 12533
a = 2;                  % approximate y-intercept for age 45 couple (p. 47)
b = 0.6;                % approximate average marg. tax rate from Table 3

rng(1);                 % set seed
norm = sort(normrnd(0,1,nAgents,1));
sim = @(th) exp(th(1) + th(2)*norm);
obj = @(th) sum((LBs - prctile(ystar(sim(th),a,b),prc)).^2);

opts = optimset('display','iter','maxit',200);
theta = ktrlink(obj,th0,[],[],[],[],[0; 0],[3; 3],[],opts)
