%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to do data analysis using sub-functions
% Outputs: Logtheta is the log of the hyperparameter of Brownian Motion, PP is the sampledd marginal firing probabilities, Samp is the sample of copula parameters, 
% logwt is the log weights of the spherical HMC samples, acpt is the acceptance rate of spherical HMC, and TIME is the time to do MCMC sampling.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Logtheta, PP, Samp, logwt, acpt, TIME] = Main( )
clear;

%% set rng seed
%s = RandStream('mt19937ar','Seed',2013);
%RandStream.setGlobalStream(s);
s = RandStream('mcg16807','Seed',2013)
RandStream.setDefaultStream(s)

%% Initialization for samplers
w = 0.5;
m = 10;
Niteration = 3000;
BurnIn=1000;

% Source Spike Train Data and Initialization
DATA = load('Data016.txt'); 
D = 9;
N = 51;
n = 2001;
time = 5*(0:2000)/1000;
x=time;
DATA1 = zeros( N, n, D );
DATA2 = zeros( N, n, D );
for( i = 1:D )    
    DATA1(:,:,i) = transpose( DATA( DATA(:,52)==i & DATA(:,53)==1, 1:51 ) ); % DATA1 is the spike train data under lever 1
    DATA2(:,:,i) = transpose( DATA( DATA(:,52)==i & DATA(:,53)==2, 1:51 ) ); % DATA2 is the spike train data under lever 2
end

Y1 = DATA1(:,:,[1:5]); % Analyzing the first 5 active neurons under lever 1
matlabpool('open',12); % Open matlab pools for parallel computing
[Logtheta, PP, Samp, logwt, acpt, TIME] = CopulaSHMC( Y1, Niteration, BurnIn, w, m );
matlabpool('close'); % Close matlab pools
save( 'Data016_Lever1_12345.mat' )

Y2 = DATA2(:,:,[1:5]); % Analyzing the first 5 active neurons under lever 2
matlabpool('open',12); % Open matlab pools for parallel computing
[Logtheta, PP, Samp, logwt, acpt, TIME] = CopulaSHMC( Y2, Niteration, BurnIn, w, m );
matlabpool('close'); % Close matlab pools
save( 'Data016_Lever2_12345.mat' )


end
