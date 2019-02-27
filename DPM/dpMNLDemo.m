% This program performs nonlinear classification using Dirichlet process mixtures. 
% We model the joint distribution of x and y as a Dirichlet process 
% mixture of some simple distributions.Within each component of the mixture, 
% we use an MNL model to capture the dependence of y on x. 

% While the program is running, predictive probabilities on the test set is
% estimated after discarding the initial samples (specified by burnIn).
% These probabilites are save in a file call "dpMnlProb.dat".

% At each iteration, the size of each component, and the probability of
% accepting new regression parameters in re-sampling will be printed on the 
% screen. To make the program run faster, you can remove all 'disp' codes
% in this program.

% After running this program, the posterior samples are saved in a .mat
% file named "dpMnlSamp.mat", and the results (i.e., the accuracy rate and 
% the F1 measure) on the test set are saved in dpMnl_results.dat.

function dpMNL_Demo(numberOfIterations, burnIn);

seedNum = 10; % This is the seed number to simulate sample data. Change this if 
% you want to create a different dataset

% This part simulates 10000 points of which only 100 will be used for training. 
% The rest are used as the test set. 
[xTrain, yTrain, yTrainD, xTest, yTest, yTestD] = simulateData(10000, seedNum); 
% xTrain: covariates for the training set   
% yTrain: response variable for the training set
% yTrainD: dummy variables for yTrain
% xTest: covariates for the test set  
% yTest: response variable for the test set
% yTestD: dummy variables for yTest

% This part sends the datasets to dpMnl function to obtain the 
% predictive probabilities on the test set.
p = dpMnl(xTrain, yTrain, yTrainD, xTest, yTest, yTestD, numberOfIterations, burnIn); 

% This part calculates the accuracy rate and the F1 measure on the test
% set.
result = getResults(yTest, p) 

% This part saves the results in a file.
dlmwrite('dpMnl_results.dat', result, '-append'); 


function predProb = dpMnl(xTrain, yTrain, yTrainD, xTest, yTest, yTestD, numberOfIterations, burnIn);

global nTrain nTest n d nLeaf nCov;
global mu0 Sigma0 mu00 Sigma00 a0 b0;
global aSigma0 bSigma0 muMu sigMu;
global aSigma00 bSigma00 muSig sigSig;
global abNuA abNuB sigBeta;

global leapFrog eps;

leapFrog = 200; % Number of steps for the Hamiltonian dynamics
eps = 0.2; % This is the constant multiplier for the step size 

[nTrain, nCov] = size(xTrain);
nTest = size(xTest, 1);
data = [xTrain; xTest];
data = data - repmat(mean(data), nTrain+nTest, 1); % This part centers the data;
xTrain = data(1:nTrain, :);
xTest  = data(nTrain+1:end, :);

[n, nCov] = size(xTrain);
nLeaf = size(yTrainD, 2);

source(1).ind = 1:nCov;
% These are the parameters of the gamma prior for scale paramtere, alpha.
a0 = -3; b0 = 2; 

% mu(j) ~ N(mu0(j), sigma0(j)), where sigma0 is the standard deviation 
% mu0(j) ~ N(muMu, sigMu) 
% log(sigma0(j).^2) ~ N(aSigma0, bSigma0)
mu0 = zeros(1, nCov);
muMu = zeros(1, nCov); sigMu = 2*ones(1, nCov);
Sigma0 = 1*ones(1, nCov);
aSigma0 = 0; bSigma0 = 1;

% log(Sigma(j)^2) ~ N(mu00(j), sigma00(j)), where Sigma00 is the standard
% deviation
% mu00(j) ~ N(muSig, sigSig)
% log(Sigma00(j).^2) ~ N(aSigma00, bSigma00)
mu00 = zeros(1, nCov);
muSig = 0; sigSig = 1;
Sigma00 = 1*ones(1, nCov);
aSigma00 = 0; bSigma00 = 1;

abNuA  = [0, 1]; % Paramters of the prior for the intercept
abNuB  = [0, 1]; % Parameters of the prior for the coefficients

% Initial values
nuA = 1;
nuB = 1;
sigBeta = ones(nCov+1, nLeaf);

sigComp = sigBeta;
sigComp(1, :) = nuA*sigBeta(1, :);
sigComp(2:end, :) = nuB*sigBeta(2:end, :);

% An initial value for the Scale parameter of the Drichlet process
alpha(1) = .001;

% thetaStar holds the unique parameters sampled from the Dirichlet process
thetaStar.mu = zeros(1, nCov);
thetaStar.sd = 2*ones(1, nCov);
thetaStar.nuA = nuA;
thetaStar.nuB = nuB;
thetaStar.sigComp = sigComp;
thetaStar.beta = ones(1+nCov, nLeaf);

% An initial value for the cluster identifiers, J
newJ = repmat(1, n, 1); 

% An initial value for the frequencey of each cluster
newNj = n;

startTime = cputime;

newProb = zeros(nTest, nLeaf);
countP = 0; 
counter = 0;
sumPX = 0;
sumProb = 0;
meanAccept = 0;
for iter = 1:numberOfIterations
iter

% This part calls the MCMC algorithm
[thetaStar, newJ, newNj] = main_MCMC(xTrain, yTrainD, thetaStar, newJ, newNj, alpha(iter));
% This prints the new frequncy of each cluster on the screen 
disp('Number of samples in each component:')
disp(newNj')
% This part resamples the parameters of the Dirichlet process given 
% the current assignment of data points to cluster.
[thetaStar, acceptP] = remix(xTrain, yTrainD, newJ, thetaStar); 
disp('Acceptance prob when resampling coefficients: ')
disp(acceptP)

% This part samples new alpha from its posterior distribution 
alpha(iter+1) = pickAlpha(alpha(iter), length(newNj), a0, b0); 

nComp = length(newNj); % Number of compontnst (i.e., clusters) in the mixture


% This parts sets the constant multiplier for the step size of the
% Hamiltonian dynamics. If the acceptance rates are very low or very high,
% you can change these values to obtain a more appropriate acceptance rate.
if rem(iter, 5) ==0 
    eps = 0.4;
else
    eps = 0.2;
end

% This part obtains the predictive probability for the test set when the
% number of iterations is larger than "burnIn" value.
if iter > burnIn
    countP = countP+1;
    [pY, pX, q] = getPredProb(xTest, thetaStar, newNj, alpha(iter+1), mu0, Sigma0, mu00, Sigma00, sigBeta);
    pX = repmat(pX, 1, nLeaf);
    sumPX = pX+sumPX;    
    p = pY.*pX;
    sumProb = (p+sumProb);    
    predProb = (sumProb./sumPX);
    dlmwrite('dpMnlProb.dat', predProb);
end


% This part samples from the posterior distribution of hyperparameters. You
% can remove this part if you do not want to use hyperparamters.
if iter >= 5

    uniqueMu = cat(1, thetaStar.mu);
    uniqueSd = log( (cat(1, thetaStar.sd)).^2 );
    
    for j = 1:nCov
        mu0(1, j) = getMu0(uniqueMu(:, j), mu0(1, j), muMu, sigMu, Sigma0(1, j));   
        Sigma0(1, j) = sqrt(getSig0(uniqueMu(:, j), Sigma0(1, j).^2, aSigma0, bSigma0, mu0(1, j)));
        
        mu00(1, j) = getMu0(uniqueSd(:, j), mu00(1, j), muSig, sigSig, Sigma00(1, j));   
        Sigma00(1, j) = sqrt(getSig0(uniqueSd(:, j), Sigma00(1, j).^2, aSigma00, bSigma00, mu00(1, j)));
    end
        
end


dpMnlSamp(iter).mu = cat(1, thetaStar.mu);
dpMnlSamp(iter).sd = cat(1, thetaStar.sd);
dpMnlSamp(iter).nuA = cat(1, thetaStar.nuA);
dpMnlSamp(iter).nuB = cat(1, thetaStar.nuB);
dpMnlSamp(iter).beta = cat(3, thetaStar.beta);

save dpMnlSamp dpMnlSamp;


end

timePassed = cputime - startTime




function [thetaStar, updatedJ, updatedNj] = main_MCMC(x, y, thetaStar, J, nj, alpha);

global a0 b0; 
global mu0 Sigma0 mu00 Sigma00; 
global aSigma00 bSigma00 muSig sigSig;
global abNuA abNuB sigBeta;
global n nTrain nTest nLeaf;
M=5; % Number of auxillary components

for i = 1:n
    curInd = J(i);        
    nj(curInd) = nj(curInd) - 1;
    if nj(curInd) == 0
        phi = thetaStar(curInd);
        nj(curInd) =[];
        thetaStar(curInd)= [];
        J(J>curInd) = J(J>curInd) - 1;
        kBar = length(nj);
        thetaStar(kBar+1) = phi;
        for m=1:(M-1)
            sd = sqrt(exp(normrnd(mu00, Sigma00)));
            mu = normrnd(mu0, Sigma0);
            nuA = sqrt(exp(normrnd(abNuA(1), abNuA(2))));
            nuB = sqrt(exp(normrnd(abNuB(1), abNuB(2))));
            sigComp = sigBeta;
            sigComp(1, :) = nuA*sigBeta(1, :);
            sigComp(2:end, :) = nuB*sigBeta(2:end, :);
            beta = normrnd(0, sigComp);            
            thetaStar(kBar+1+m).mu = mu;
            thetaStar(kBar+1+m).sd = sd;
            thetaStar(kBar+1+m).nuA = nuA;
            thetaStar(kBar+1+m).nuB = nuB;
            thetaStar(kBar+1+m).sigComp = sigComp;
            thetaStar(kBar+1+m).beta = beta;
        end
    else
        kBar = length(nj);
        for m=1:M
            sd = sqrt(exp(normrnd(mu00, Sigma00)));
            mu = normrnd(mu0, Sigma0);
            nuA = sqrt(exp(normrnd(abNuA(1), abNuA(2))));
            nuB = sqrt(exp(normrnd(abNuB(1), abNuB(2))));
            sigComp = sigBeta;
            sigComp(1, :) = nuA*sigBeta(1, :);
            sigComp(2:end, :) = nuB*sigBeta(2:end, :);
            beta = normrnd(0, sigComp);            
            thetaStar(kBar+m).mu = mu;
            thetaStar(kBar+m).sd = sd;
            thetaStar(kBar+m).nuA = nuA;
            thetaStar(kBar+m).nuB = nuB;
            thetaStar(kBar+m).sigComp = sigComp;
            thetaStar(kBar+m).beta = beta;
         end
    end

    q1 = zeros(kBar, 1);
        
    for k = 1:kBar
        mu = thetaStar(k).mu;
        beta = thetaStar(k).beta;
        eta = [1, x(i, :)]*beta;
        sd = thetaStar(k).sd;
        sigComp = thetaStar(k).sigComp;
        q1(k) = getLogLike(x(i, :), y(i, :), mu, sd, eta);
    end
    q1 = q1 + log(nj) - log(n-1+alpha);
    
    q2 = zeros(M, 1);
    for k = 1:M
        mu = thetaStar(kBar+k).mu;
        beta = thetaStar(kBar+k).beta;
        eta = [1, x(i, :)]*beta;
        sd = thetaStar(kBar+k).sd;
        sigComp = thetaStar(kBar+k).sigComp;
        q2(k) = getLogLike(x(i, :), y(i, :), mu, sd, eta);
    end

    q2 = q2+(log(alpha) - log(M))-log(n-1+alpha);
    q = [q1; q2];

    qMax = max(q);
    qRel = q - qMax;
    
    q = exp(qRel);
    
    q = q./sum(q);

    qCumSum = repmat(0, length(q), 1);
    qCumSum = cumsum (q); 

    u = rand;
    k0 = find ( qCumSum >= u);    
    picked = k0(1);
    
    if picked <= kBar
        J(i) = picked;
        nj(picked) = nj(picked)+1;
        thetaStar(kBar+1:end) = [];        
    else
        J(i) = kBar+1;
        nj = [nj; 1];
        phi = thetaStar(picked);
        thetaStar(kBar+1:end) = [];
        thetaStar(kBar+1) = phi;
    end

end    
    
updatedNj = nj;
updatedJ = J;



function result = getLogLike(x, y, mu, sd, eta);
global n d nLeaf;

sd = sd + .001;
diffX = (x - mu);
logLikeX = sum( - log(sd) - ( ((diffX).^2)./(2*(sd.^2)) ) ); 

m = max(eta');
logLikeY = sum( (sum((y.*eta), 2) - (m' + log(sum(exp(eta - repmat(m', 1, nLeaf) ), 2)) ) ) );

result = logLikeX + logLikeY;


function [remixedThetaStar, acceptProb] = remix(x, y, J, thetaStar)
 
global mu0 Sigma0 mu00 Sigma00 abNuA abNuB sigBeta n d nLeaf nCov;
global eps;

u1 = sort(unique(J));
for i = 1:length(u1)
    X = x(J == i, :);
    Y = y(J == i, :);
    nY = length(J(J==i));

    mu = thetaStar(i).mu;
    sd = thetaStar(i).sd;
    nuA = thetaStar(i).nuA;
    nuB = thetaStar(i).nuB;
    sigComp = thetaStar(i).sigComp;
    beta=thetaStar(i).beta;

    for j = 1:nCov
        mu(j) = getMu0(X(:, j), mu(j), mu0(j), Sigma0(j), sd(j));   
        sd(j) = sqrt(getSig0(X(:, j), sd(j).^2, mu00(j), Sigma00(j), mu(j)));
    end

    relatedBeta = reshape(beta(1, :), 1, nLeaf);
    nuA = sqrt(getSigBeta(relatedBeta, nuA^2, abNuA(1), abNuA(2)));
    relatedBeta = reshape(beta(2:end, :), 1, nLeaf*nCov);
    for repHype = 1:50
    nuB = sqrt(getSigBeta(relatedBeta, nuB^2, abNuB(1), abNuB(2)));
    end
    sigComp = sigBeta; 
    sigComp(1, :) = nuA*sigComp(1, :);
    sigComp(2:end, :) = nuB*sigComp(2:end, :);


    e = eps*( 1 ./ sqrt( 1 ./ (sigComp .^ 2) + nY/4) );

    [beta, acceptProb(1, i)] = getBeta([ones(nY, 1), X], Y, beta, sigComp, e);
    
    
    remixedThetaStar(i).mu = mu;
    remixedThetaStar(i).beta = beta;
    remixedThetaStar(i).sd = sd;
    remixedThetaStar(i).nuA = nuA;
    remixedThetaStar(i).nuB = nuB;
    remixedThetaStar(i).sigComp = sigComp;
end



function [predProb, postX, q] = getPredProb(x, thetaStar, nj, alpha, mu0, Sigma0, mu00, Sigma00, sigBeta);

global a0 b0; 
global aSigma00 bSigma00 muSig sigSig abNuA abNuB;

global n nTrain nTest nLeaf nCov;

iStar = length(nj);
q = zeros(nTest, iStar+1);

% This part is sample from G_0, which is used to get the predictive
% probability
sd0 = repmat(sqrt(exp(normrnd(mu00, Sigma00))), nTest, 1);
mu0 = repmat(normrnd(mu0, Sigma0), nTest, 1);
nuA0 = sqrt(exp(normrnd(abNuA(1), abNuA(2))));
nuB0 = sqrt(exp(normrnd(abNuB(1), abNuB(2))));
sigComp0 = sigBeta;
sigComp0(1, :) = nuA0*sigBeta(1, :);
sigComp0(2:end, :) = nuB0*sigBeta(2:end, :);
beta0 = normrnd(0, sigComp0);
thetaStar(iStar+1).beta = beta0;

% This part uses the unique component to get the predictive probability
for k = 1:iStar
    mu = repmat(thetaStar(k).mu, nTest, 1);
    sd = repmat(thetaStar(k).sd, nTest, 1);
    diffX = (x - mu);
    q(:, k) = sum( - log(sd) - ( ((diffX).^2)./(2*(sd.^2)) ), 2);     
    q(:, k) = q(:, k) + log(nj(k)) - log(n+alpha);    
end

diffX = (x - mu0);
q(:, iStar+1) = sum( - log(sd0) - ( ((diffX).^2)./(2*(sd0.^2)) ), 2);     
q(:, iStar+1) = q(:, iStar+1) + log(alpha) - log(n+alpha);

% This calculates P(x) for the current iteration
m = max(q');
postX(:, 1) = exp(m'+log(sum(exp(q - repmat(m', 1, iStar+1) ), 2)));


% This part gets P(y) for the current iteration
qRel = q - repmat(m', 1, iStar+1);
q = exp(qRel);
q = q./repmat(sum(q, 2), 1, iStar+1);

for k = 1:iStar+1
    etaT = [ones(nTest, 1), x]*thetaStar(k).beta;
    m = max(etaT');        
    predProb(:, :, k) = repmat(q(:, k), 1, nLeaf).*exp(etaT - repmat(m' + log(sum(exp(etaT - repmat(m', 1, nLeaf) ), 2)), [1, nLeaf, 1] ) );
end
    
predProb = sum(predProb, 3);




%******************* PickAlpha ******************************

% This function updates "alpha: based on the old alpha and number of unique
% theta's

function newAlpha = pickAlpha(oldAlpha, iStar, a0, b0)
m = 40; w = 5;

x = log(oldAlpha+0.001); 
z = getLogPostAlpha(x, iStar, a0, b0) - exprnd(1);

u = rand;
L = x - w*u;
R = L + w;
v = rand;
J = floor(m*v);
K = (m-1) - J;

while J>0 && z < getLogPostAlpha(L, iStar, a0, b0)
    L = L - w;
    J = J - 1;
end

while K>0 && z < getLogPostAlpha(R, iStar, a0, b0)
    R = R+w;
    K = K-1;
end

u = rand;
newX = L + u*(R-L);

while z > getLogPostAlpha(newX, iStar, a0, b0)
    if newX < x
        L = newX;
    else
        R = newX;
    end    
    u = rand();
    newX = L + u*(R-L);
end

newAlpha = exp(newX);

    
function logPost = getLogPostAlpha(x, iStar, a0, b0)
global n;
alpha = (exp(x));

logLike = iStar*log(alpha) + gammaln(alpha) - gammaln(alpha+n); 

logPrior = sum(-.5*log(2*pi)-log(b0)-( (( x - a0 ).^2)./(2*(b0.^2)) ));

logPost = logLike + logPrior;


function mewAlpha = pickAlpha0 (oldAlpha, iStar)


global a0 b0 n;

nu = betarnd(oldAlpha+1, n);

proportion = (a0 + iStar - 1)/(n*(b0-log(nu)));
piNu = proportion /(proportion +1);

u = rand;
if u <= piNu
    mewAlpha = gamrnd(a0+iStar, 1./(b0-log(nu)));
else
    mewAlpha = gamrnd(a0+iStar-1, 1./(b0-log(nu)));
end



function newMu = getMu0(y, oldMu, mu0, sigma0, sigma)
m = 40; w = 5;

x = oldMu; 
z = getLogPostMu0(y, x, mu0, sigma0, sigma) - exprnd(1);

u = rand();
L = x - w*u;
R = L + w;
v = rand();
J = floor(m*v);
K = (m-1) - J;

while J>0 && z < getLogPostMu0(y, L, mu0, sigma0, sigma)
    L = L - w;
    J = J - 1;
end

while K>0 && z < getLogPostMu0(y, R, mu0, sigma0, sigma) 
    R = R+w;
    K = K-1;
end

u = rand();
newX = L + u*(R-L);

while z > getLogPostMu0(y, newX, mu0, sigma0, sigma)
    if newX < x
        L = newX;
    else
        R = newX;
    end
    
    u = rand();
    newX = L + u*(R-L);
end

newMu = newX;

    
function logPost = getLogPostMu0(y, mu, mu0, Sigma0, sd)

sd = sd+0.001;
Sigma0 = Sigma0+0.001;
diffY = (y - mu);
logLike = sum( -.5*log(2*pi) - log(sd) - ( ((diffY).^2)./(2*(sd.^2)) ) ); 

diffMu = mu - mu0;
logPriorMu = sum(-.5*log(2*pi)-log(Sigma0)-( (( diffMu ).^2)./(2*(Sigma0.^2)) ) );

logPost = logLike + logPriorMu; 




function newSigma2 = getSig0(y, oldSigma2, mu00, sigma00, muY)
m = 40; w = 5;

x = log(oldSigma2+0.001); 
z = getLogPostSigma0(y, x, mu00, sigma00, muY) - exprnd(1);

u = rand();
L = x - w*u;
R = L + w;
v = rand();
J = floor(m*v);
K = (m-1) - J;

while J>0 && z < getLogPostSigma0(y, L, mu00, sigma00, muY)
    L = L - w;
    J = J - 1;
end

while K>0 && z < getLogPostSigma0(y, R, mu00, sigma00, muY)
    R = R+w;
    K = K-1;
end

u = rand();
newX = L + u*(R-L);

while z > getLogPostSigma0(y, newX, mu00, sigma00, muY)
    if newX < x
        L = newX;
    else
        R = newX;
    end
    
    u = rand();
    newX = L + u*(R-L);
end

newSigma2 = exp(newX);

    
function logPost = getLogPostSigma0(y, x, mu00, Sigma00, mu)

sd = sqrt(exp(x));
sd = sd+0.001;
Sigma00 = Sigma00+0.001;

diffY = (y - mu);
logLike = sum( -.5*log(2*pi) - log(sd) - ( ((diffY).^2)./(2*(sd.^2)) ) ); 

diffSig = x - mu00;
logPriorSd = sum(-.5*log(2*pi)-log(Sigma00 )-( (( diffSig ).^2)./(2*(Sigma00.^2)) ));

logPost = logLike + logPriorSd;



function [updateB, acceptProb] = getBeta(d, r, b, sigma0, e)
global leapFrog;

[dim1, dim2] = size(b);

oldB = b;

    E = getE(d, r, oldB, sigma0);
    g = getG(d, r, oldB, sigma0);
    p = randn(size(b));
    H = sum(sum( .5*(p.^2) ))+ E;
    newB = b; 
    newG = g; 

    for leap=1:leapFrog
        p = p - e.*newG/2;
        newB = newB + e.*p;
        tempG = getG(d, r, newB, sigma0);
        if isfinite(tempG)
            newG = tempG;
        end

        p = p - e.*newG/2;
    end

    newE = getE(d, r, newB, sigma0);
    newH = sum(sum( .5*(p.^2) )) + newE;

    acceptProb = min(1, exp(H - newH)); 

    if (rand < acceptProb)
        updateB = newB;
    else
        updateB = oldB;
    end


% This calculates the energy function for the posterior distribution. 

function E = getE(d, r, beta, sigma)
global nLeaf;

eta = d*beta;

m = max(eta');
logLike = sum( (sum((r.*eta), 2) - (m' + log(sum(exp(eta - repmat(m', 1, nLeaf) ), 2)) ) ) );

logPrior =  sum(sum( ( -0.5*(beta).^2 ) ./ (sigma.^2) ));

E = -(logLike + logPrior);



% This part calculates the derivatives for all parameters. 
% Note that the calculation is completely vectorized and quite fast. 
% Moreover, the code is written in a way that avoids overflow.

function g = getG(d, r, beta, sigma)
global nLeaf;
eta = d*beta;

m = max(eta');

dLogLike = (d'*r - d'*exp( eta - repmat( (m' + log(sum(exp(eta - repmat(m', 1, nLeaf) ), 2)) ) , 1, nLeaf) ) );
dLogPrior =  -(beta) ./ (sigma.^2);
 
g = -(dLogLike + dLogPrior);





function newSigma = getSigBeta(beta, oldSigma2, mu, tau)
m = 40; w = 5;

x = log(oldSigma2); 
z = getLogPostSigma(beta, x, mu, tau) - exprnd(1);

u = rand();
L = x - w*u;
R = L + w;
v = rand();
J = floor(m*v);
K = (m-1) - J;

while J>0 && z < getLogPostSigma(beta, L, mu, tau)
    L = L - w;
    J = J - 1;
end

while K>0 && z < getLogPostSigma(beta, R, mu, tau) 
    R = R+w;
    K = K-1;
end

u = rand();
newX = L + u*(R-L);

while z > getLogPostSigma(beta, newX, mu, tau)
    if newX < x
        L = newX;
    else
        R = newX;
    end
    
    u = rand();
    newX = L + u*(R-L);
end

newSigma = exp(newX);


function logPost = getLogPostSigma(beta, x, mu, tau)

sigma = sqrt(exp(x));
n = length(beta);

logPost = -n*(log(sqrt(2*pi)) + log(sigma)) - sum( (beta).^2 )/(2*(sigma^2)) - (log(sqrt(2*pi)) + log(tau)) - ( (x - mu)^2 )/(2*(tau^2));






function results = getResults(target, p);

n = length(target);
c = unique(target);

for i = 1:n
    classProb(i) = p(i, target(i));
end
avgLogProb = mean(log(classProb));

[maxPred, maxInd] = max(p'); 
predClass = maxInd';

accRateTest = mean(logical(predClass == target))

[m, indM] = max(p');
predClass = indM';


for j = 1:length(c)
    i = c(j);
    categA = length(find(target==i & predClass == i));
    categB = length(find(target~=i & predClass == i));
    categC = length(find(target==i & predClass ~= i));
    
    categF1(i) = 2*categA / (2*categA+ categB + categC);
end

f1 = mean(categF1);


for i = 1:n
    precisionI(i) = 1 / sum(logical(p(i, :) >= p(i, target(i))));
end

precision = mean(precisionI);


results = [avgLogProb, accRateTest, precision, f1];   




function [train, rTrain, rTrainD, test, rTest, rTestD] = simulateData(n, state)

rand('state', state);
randn('state', state);

nClass = 4;
nVar = 5;

mu = normrnd(0, 1, 2, nVar);
sd = sqrt(exp(normrnd(0, 2, 2, nVar)));

sigmaInt = .1;
nu = sqrt(exp(normrnd(0, 2, 1, 2)));

beta1 = [normrnd(0, sigmaInt, 1, nClass); normrnd(0, nu(1), nVar, nClass)];
beta2 = [normrnd(0, sigmaInt, 1, nClass); normrnd(0, nu(2), nVar, nClass)];

x1 = normrnd(repmat(mu(1, :), n/2, 1), repmat(sd(1, :), n/2, 1), n/2, nVar);
etaT = [ones(n/2, 1), x1]*beta1;
etaT = etaT + normrnd(0, .05, size(etaT));
p = exp(etaT - repmat(log(sum(exp(etaT), 2)), [1, nClass, 1]));   
[maxPred, maxInd] = max(p'); 
y1 = maxInd';

x2 = normrnd(repmat(mu(2, :), n/2, 1), repmat(sd(2, :), n/2, 1), n/2, nVar);
etaT = [ones(n/2, 1), x2]*beta2;
etaT = etaT + normrnd(0, .05, size(etaT));
p = exp(etaT - repmat(log(sum(exp(etaT), 2)), [1, nClass, 1]));   
[maxPred, maxInd] = max(p'); 
y2 = maxInd';

x = [x1; x2];
y = [y1; y2];

yD = zeros(n, nClass);
for i = 1:n
    yD(i, y(i)) = 1;
end

samp = randsample(n, 100);
train = x(samp, :);
rTrain = y(samp, :);
rTrainD = yD(samp, :);
test = x(setdiff(1:n, samp), :);
rTest = y(setdiff(1:n, samp), :);
rTestD = yD(setdiff(1:n, samp), :);

% gscatter(train(:, 2), train(:, 3), rTrain, '', '+*>d');