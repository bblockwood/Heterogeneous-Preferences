# AMPL model file


################################################################################
##  Define parameters              

# Utility parameters
param sigma;
param gamma;

# Characterize agents
param nAgents;
set I := 1..nAgents;
param mu;
param sd;
param phibar;
param yLF {I} := exp(Normal(mu,sd));
param lambda {i in I} := yLF[i]^sigma;
param phi {I} := phibar;
param theta {i in I} := lambda[i]^(phi[i]); 
param w {i in I} := lambda[i]^((1-phi[i])/sigma);

# Characterize planner: preference neutral cardinalization
param alphaWeights {i in I} = (theta[i])^((1-sigma)/(sigma+gamma-1));


################################################################################
##  Define choice variables              

# Individual chooses earnings (initial guess = laissez faire earnings)
var y {i in I} >= 0 := yLF[i];

# Govt chooses linear tax regime: c = a + b*y
var a;
var b >= 0;

# Earnings and tax regime pin down c, l, and therefore u
var c {i in I} = a + b*y[i];
var l {i in I} = y[i]/w[i];
var u {i in I} = if gamma = 1 
	then theta[i]*log(c[i]) - l[i]^sigma/sigma  # log utility if gamma == 1
	else theta[i]*(c[i]^(1-gamma)-1)/(1-gamma) - l[i]^sigma/sigma;


################################################################################
##  Define optimization problem: maximize swf s.t. agents' FOCs & govt BC       

maximize swf: sum {i in I} alphaWeights[i]*u[i];

subject to 
agentFOC {i in I}: b*lambda[i] = c[i]^gamma*y[i]^(sigma-1);
govtBC: sum {i in I} (-a + (1-b)*y[i]) = 0;

problem SwfMaximization:
	swf,
	y, a, b, c, l, u,
	agentFOC, govtBC;