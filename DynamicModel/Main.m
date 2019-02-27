
clear;
s = RandStream('mcg16807','Seed',2013);
RandStream.setDefaultStream(s)


%% Initialization for samplers
w = 0.5;
m = 10;
Niteration = 3000;
BurnIn=1000;

%% Real Data Analysis
DATA = load('Data016_100ms.txt');  %% Load dataset from text files
D = 9; %% Number of neurons
N = 51; %% Number of trials
n = 101; %% Number of time bins
time = (0:100)/10;
DATA1 = zeros( N, n, D ); 
DATA2 = zeros( N, n, D ); 
Y1 = zeros( D, n, N); %% Data for scenario 1: rewarded stimulus
Y2 = zeros( D, n, N); %% Data for scenario 2: non-rewarded stimulus
Y1temp = zeros( N, n );
Y2temp = zeros( N, n );
for( i = 1:D )    
    DATA1(:,:,i) = transpose( DATA( DATA(:,52)==i & DATA(:,53)==1, 1:51 ) );
    DATA2(:,:,i) = transpose( DATA( DATA(:,52)==i & DATA(:,53)==2, 1:51 ) );
    Y1temp(:,:) = DATA1(:,:,i);
    Y1(i,:,:) = transpose(Y1temp);
    Y2temp(:,:) = DATA2(:,:,i);
    Y2(i,:,:) = transpose(Y2temp);   
end



matlabpool('open',2);  
[Logeta, Logrho, Logalpha, LogJ2, TR, UU, II, SIGMA, TIME] = MCMC( time, Y1([1,7],:,:), Niteration, w, m );
matlabpool('close')

save( 'Rslt.mat' )

