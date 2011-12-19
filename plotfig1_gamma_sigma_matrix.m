% Plots output from tastes.run AMPL file. 

clear all;
clc;

% Load data
dataFile = importdata('mtr.csv');
assignin('base', 'gamma', dataFile.data(:,1));
assignin('base', 'sigma', dataFile.data(:,2));
assignin('base', 'phi', dataFile.data(:,3));
assignin('base', 'mtr', dataFile.data(:,4));

gammas = unique(gamma);
sigmas = unique(sigma);
nGammas = size(gammas,1);
nSigmas = size(sigmas,1);

yMin = -0.5;
yMax = 1;
axisBnds = [min(phi)
            max(phi)
            yMin
            yMax];

i = 1; % index
for g=1:nGammas
    for s=1:nSigmas
        % Restrict to phis and mtrs for this gamma and sigma
        phiArray = phi(gamma==gammas(g) & sigma==sigmas(s));
        mtrArray = mtr(gamma==gammas(g) & sigma==sigmas(s));
        
        subplot(nGammas,nSigmas,i);
        hold on;
        line([0 0],[yMin yMax]);
        line([1 1],[yMin yMax]);
        line(phiArray,mtrArray,'Color',[0 0 1],'LineWidth',1.5);%figure(gcf)
        axis(axisBnds);
        hold off;

        if mod(i,3)==1, ylabel('Marginal tax rate'), end
        if i>(nGammas-1)*nSigmas, xlabel('\phi'), end
        title(['\gamma=' num2str(gammas(g)) ', ' '\sigma=' num2str(sigmas(s))]);

        
        i = i+1;
    end
end