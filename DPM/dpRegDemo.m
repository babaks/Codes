% The only inputs required are numIter (number of Iterations)
% burnIn (number of initial MCMC samples that should be ignored). For this
% problem, you can set numIter = 500; burnIn = 200.

function dpRegDemo(numIter, burnIn)

% This part generates the data
n = 1050; nTrain = 50; i = 1;
[xTrain, yTrain, xTest, yTest] = createData(n, nTrain, i); 


% This part calles the main program. The outputs are the number of clusters
% (i.e., components of the mixture) and posterior samples of parameters.
[samp, nClus, ASE] = dpLR(xTrain, yTrain, xTest, yTest, numIter, burnIn);
disp('The average squared error on the test set is:')
disp(ASE)

J = cat(2, samp.J);
for i = 1:nTrain
    [f, c] = hist(J(i, :), 1:nClus);
    dpProb(i, :) = f/sum(f);
end
[maxP, dpid] = max(dpProb, [], 2);

% This part gets the posterior expectation of regression parameters for
% each component.
for i = 1:length(samp)
theta(i).beta = cat(2, samp(i).thetaStar(1:nClus).beta);
end
meanBeta = mean(cat(3, theta.beta), 3);


% This part plots the final results
if nClus > 1
    gscatter(xTrain, yTrain, dpid)
    hold on
    for j = 1:nClus
    xj = xTrain(dpid==j);    
    x0 = [min(xj):.2:max(xj)];
    yHat = [ones(length(x0), 1), x0']*meanBeta(:, j);
    plot(x0, yHat, '-g');
    end
else
    scatter(xTrain, yTrain)
    hold on    
    x0 = [min(x):.2:max(x)];
    yHat = [ones(length(x0), 1), x0']*meanBeta;
    plot(x0, yHat, '-g');
end


% This is the main program
function [samp, nClus, meanSqErr] = dpLR(xTrain, yTrain, xTest, yTest, numberOfIterations, burnIn)
global nTrain nTest n nCov nCont nBi nPo;
global mu0 Sigma0 mu00 Sigma00 aBi bBi aPo bPo a0 b0;
global aSigma0 bSigma0 muMu sigMu;
global aSigma00 bSigma00 muSig sigSig;
global aSigmaY bSigmaY;
global abNuA abNuB sigBeta;
global misFlagTrain misFlagTest contFlag poissFlag binFlag;
global leapFrog epsVec;

leapFrog = 100; % Number of steps for the Hamiltonian dynamics
epsVec = [0.1, 0.3];  

[nTrain, nCov] = size(xTrain);
nTest = size(xTest, 1);


n = nTrain;

source(1).ind = 1:nCov;


contFlag = [1:nCov];
poissFlag = [];
binFlag = [];

nCont = length(contFlag);
nPo = length(poissFlag);
nBi = length(binFlag);


misFlagTrain = ones(size(xTrain));
misFlagTrain(isnan(xTrain)) = 0; 

xTrain(isnan(xTrain)) = 0;

misFlagTest = ones(size(xTest));
misFlagTest(isnan(xTest)) = 0; 

xTest(isnan(xTest)) = 0;

% These are the parameters of the gamma prior for scale paramtere, alpha.
a0 = -3; b0 = 2; 

% mu(j) ~ N(mu0(j), sigma0(j)), where sigma0 is the standard deviation 
% mu0(j) ~ N(muMu, sigMu) 
% log(sigma0(j).^2) ~ N(aSigma0, bSigma0)
mu0 = zeros(1, nCont);
muMu = zeros(1, nCont); sigMu = 3*ones(1, nCont);
Sigma0 = 3*ones(1, nCont);
aSigma0 = 0; bSigma0 = 1;

% log(Sigma(j)^2) ~ N(mu00(j), sigma00(j)), where Sigma00 is the standard
% deviation
% mu00(j) ~ N(muSig, sigSig)
% log(Sigma00(j).^2) ~ N(aSigma00, bSigma00)
mu00 = zeros(1, nCont);
muSig = 0; sigSig = 1;
Sigma00 = 1*ones(1, nCont);
aSigma00 = 0; bSigma00 = 1;

% muPo ~ gamma(aPo, bPo)
aPo = .5*ones(1, nPo); bPo = 10*ones(1, nPo);

% muBi ~ beta(aBi, bBi) 
aBi = 1*ones(1, nBi); bBi = 1*ones(1, nBi);

SigmaY = 1;
aSigmaY = 0; bSigmaY = 2;

abNuA  = [0, 1]; % Paramters of the prior for the intercept
abNuB  = [0, 1]; % Parameters of the prior for the coefficients

% Initial values
nuA = 3;
nuB = 3;
sigBeta = ones(nCov+1, 1);

sigComp = sigBeta;
sigComp(1, :) = nuA*sigBeta(1, :);
sigComp(2:end, :) = nuB*sigBeta(2:end, :);

% An initial value for the Scale parameter of the Drichlet process
alpha(1) = .001;

% thetaStar holds the unique parameters sampled from the Dirichlet process
thetaStar.mu = zeros(1, nCont);
thetaStar.sd = 2*ones(1, nCont);
thetaStar.sdy = 2;
thetaStar.nuA = nuA;
thetaStar.nuB = nuB;
thetaStar.sigComp = sigComp;
thetaStar.beta = ones(1+nCov, 1);
thetaStar.muPo = 1*ones(1, nPo);
thetaStar.muBi = .5*ones(1, nBi);


% An initial value for the cluster identifiers, J
newJ = repmat(1, n, 1); 

% An initial value for the frequencey of each cluster
newNj = n;

startTime = cputime;

countP = 0; 
countS = 0;
sumPX = 0;
sumProb = 0;
% 
% fig=figure;
% set(fig,'DoubleBuffer','on');
% set(gca, 'xlim',[min(xTest)-0.5 max(xTest)+0.5], 'ylim',[min(yTest)-.5 max(yTest)+.5],...
%     'NextPlot','replace','Visible','off');
% mov = avifile('dpRegMovie0.avi', 'fps', 4);
% 
% plot(xTest, yTest, '.');
% axis tight
% 
% xlabel('x');
% ylabel('y');
%         
% rect = [0 0 435 343];

% % F = getframe(gca, rect);
% % mov = addframe(mov,F);

pause(1);

for iter = 1:numberOfIterations
iter

% % This part calls the MCMC algorithm
[thetaStar, newJ, newNj] = main_MCMC(xTrain, yTrain, thetaStar, newJ, newNj, alpha(iter));


% % This prints the new frequncy of each cluster on the screen 
% disp('Number of samples in each component:')
% disp(newNj')
% % This part resamples the parameters of the Dirichlet process given 
% % the current assignment of data points to cluster.

thetaStar = remix(xTrain, yTrain, newJ, thetaStar); 

% disp('Acceptance prob when resampling coefficients: ')
% disp(acceptP)

% This part samples new alpha from its posterior distribution 
alpha(iter+1) = pickAlpha(alpha(iter), length(newNj), a0, b0); 
newNj'
nComp = length(newNj); % Number of compontnst (i.e., clusters) in the mixture


% This part obtains the predictive probability for the test set when the
% number of iterations is larger than "burnIn" value.
if iter > burnIn
    [pY, pX, q] = getPredProb(xTest, thetaStar, newNj, alpha(iter+1));
%    pX = repmat(pX, 1, nLeaf);
    sumPX = pX+sumPX;    
    
    p = pY.*pX;

    if ~isnan(sum(sum(p)))
    countP = countP+1;
    sumProb = (p+sumProb);    
    predY = (sumProb./sumPX);
    end

%     dlmwrite('dpRegPred.dat', predY);
%     mse(countP, 1) = mean((yTest - predY).^2);
end


% This part samples from the posterior distribution of hyperparameters. You
% can remove this part if you do not want to use hyperparamters.
if iter >= 5

    uniqueMu = cat(1, thetaStar(newJ).mu);
    uniqueSd = log( (cat(1, thetaStar(newJ).sd)).^2 );
    
    for j = 1:nCont
        mu0(1, j) = getMu0(uniqueMu(:, j), mu0(1, j), muMu, sigMu, Sigma0(1, j));   
        Sigma0(1, j) = sqrt(getSig0(uniqueMu(:, j), Sigma0(1, j).^2, aSigma0, bSigma0, mu0(1, j)));
        
%         mu00(1, j) = getMu0(uniqueSd(:, j), mu00(1, j), muSig, sigSig, Sigma00(1, j));   
%         Sigma00(1, j) = sqrt(getSig0(uniqueSd(:, j), Sigma00(1, j).^2, aSigma00, bSigma00, mu00(1, j)));
    end
        
        
end


if iter >= burnIn
    countS = countS + 1;
    samp(countS).thetaStar = thetaStar;
    samp(countS).newNj = newNj;
    samp(countS).J = newJ;    
    samp(countS).alpha = alpha(iter+1);
    sampComp(countS) = nComp;
end


    if length(newNj)==1

        beta = cat(2, thetaStar.beta);
        
        plot(xTrain, yTrain, 's')
        axis tight        
        hold on
        x0 = [min(xTrain):.2:max(xTrain)]';
        yHat0 = [ones(length(x0), 1), x0]*beta(:, 1);
        plot(x0, yHat0, '-g')
        
%         if iter > 50
%         F = getframe(gca, rect);
%         mov = addframe(mov,F);
%         end
        hold off 

        drawnow;
        

    elseif length(newNj)==2
        beta = cat(2, thetaStar.beta);
        ind1 = find(newJ == 1);
        ind2 = find(newJ == 2);
        plot(xTrain(ind1, :), yTrain(ind1, :), 's')
        hold on
        plot(xTrain(ind2, :), yTrain(ind2, :), 'or')               
        axis tight
        x0 = [min(xTrain(ind1, 1)):.2:max(xTrain(ind1, :))]';
        yHat0 = [ones(length(x0), 1), x0]*beta(:, 1);
        plot(x0, yHat0, '-g')

        x0 = [min(xTrain(ind2, 1)):.2:max(xTrain(ind2, :))]';
        yHat0 = [ones(length(x0), 1), x0]*beta(:, 2);
        plot(x0, yHat0, '-g')
        checkThis = 1;

        if iter > burnIn+1
            temp = sortrows([xTest, predY]);
            plot(temp(:, 1), temp(:, 2), '--g')
        end

%         if iter > 50
%         F = getframe(gca, rect);
%         mov = addframe(mov,F);
%         end
        hold off 

        drawnow;
    end


end

% mov = close ( mov );

hold off
temp = sortrows([xTest, predY]);
plot(temp(:, 1), temp(:, 2), '--g')
hold on

nClus = mode(sampComp)

timePassed = cputime - startTime

meanSqErr = mean((yTest - predY).^2);



function [thetaStar, updatedJ, updatedNj] = main_MCMC(x, y, thetaStar, J, nj, alpha)

global mu0 Sigma0 mu00 Sigma00 aPo bPo aBi bBi aSigmaY bSigmaY;
global abNuA abNuB sigBeta;
global n;
global misFlagTrain;

M=10; % Number of auxillary components

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
            muPo = gamrnd(aPo, bPo);
            muBi = betarnd(aBi, bBi);            
            sdy = sqrt(exp(normrnd(aSigmaY, bSigmaY)));
            nuA = sqrt(exp(normrnd(abNuA(1), abNuA(2))));
            nuB = sqrt(exp(normrnd(abNuB(1), abNuB(2))));
            sigComp = sigBeta;
            sigComp(1, :) = nuA*sigBeta(1, :);
            sigComp(2:end, :) = nuB*sigBeta(2:end, :);
            beta = normrnd(0, sigComp);            
            thetaStar(kBar+1+m).mu = mu;
            thetaStar(kBar+1+m).sd = sd;
            thetaStar(kBar+1+m).muPo = muPo;
            thetaStar(kBar+1+m).muBi = muBi;            
            thetaStar(kBar+1+m).sdy = sdy;
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
            muPo = gamrnd(aPo, bPo);
            muBi = betarnd(aBi, bBi);                        
            sdy = sqrt(exp(normrnd(aSigmaY, bSigmaY)));
            nuA = sqrt(exp(normrnd(abNuA(1), abNuA(2))));
            nuB = sqrt(exp(normrnd(abNuB(1), abNuB(2))));
            sigComp = sigBeta;
            sigComp(1, :) = nuA*sigBeta(1, :);
            sigComp(2:end, :) = nuB*sigBeta(2:end, :);
            beta = normrnd(0, sigComp);            
            thetaStar(kBar+m).mu = mu;
            thetaStar(kBar+m).sd = sd;
            thetaStar(kBar+m).muPo = muPo;
            thetaStar(kBar+m).muBi = muBi;                        
            thetaStar(kBar+m).sdy = sdy;
            thetaStar(kBar+m).nuA = nuA;
            thetaStar(kBar+m).nuB = nuB;
            thetaStar(kBar+m).sigComp = sigComp;
            thetaStar(kBar+m).beta = beta;
         end
    end

    q1 = zeros(kBar, 1);
        
    for k = 1:kBar
        mu = thetaStar(k).mu;
        muPo = thetaStar(k).muPo;
        muBi = thetaStar(k).muBi;
        sd = thetaStar(k).sd;
        sdy = thetaStar(k).sdy;        
        beta = thetaStar(k).beta;
        eta = [1, x(i, :)]*beta;
        q1(k) = getLogLike(x(i, :), y(i, :), mu, sd, muPo, muBi, sdy, eta, misFlagTrain(i, :));
    end
    q1 = q1 + log(nj) - log(n-1+alpha);
    
    q2 = zeros(M, 1);
    for k = 1:M
        mu = thetaStar(kBar+k).mu;
        muPo = thetaStar(kBar+k).muPo;
        muBi = thetaStar(kBar+k).muBi;
        sd = thetaStar(kBar+k).sd;
        sdy = thetaStar(kBar+k).sdy;        
        beta = thetaStar(kBar+k).beta;
        eta = [1, x(i, :)]*beta;
        q2(k) = getLogLike(x(i, :), y(i, :), mu, sd, muPo, muBi, sdy, eta, misFlagTrain(i, :));
    end

    q2 = q2+(log(alpha) - log(M))-log(n-1+alpha);
    q = [q1; q2];

    qMax = max(q);
    qRel = q - qMax;
    
    q = exp(qRel);
    
    if isnan(sum(q))
        check = 1;
    end
    q(isnan(q)) = 0;
    
    q = q./sum(q);

    qCumSum = repmat(0, length(q), 1);
    qCumSum = cumsum (q); 

    u = rand;
    k0 = find ( qCumSum >= u);    

    if numel(k0)==0
        picked = 1;        
    else
        picked = k0(1);
    end
    
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



function result = getLogLike(x, y, mu, sd, muPo, muBi, sdy, eta, misFlag)
global contFlag poissFlag binFlag;

xCont = x(:, contFlag);
xPo = x(:, poissFlag);
xBi = x(:, binFlag);

sd = sd + .001;
sdy = sdy + .001;
diffX = (xCont - mu);

logLikeXcont = sum(sum( ( - log(sd) - ( ((diffX).^2)./(2*(sd.^2)) ) ) .* misFlag(:, contFlag) )); 

logLikeXpo = sum(sum( ( log(muPo).*xPo - muPo ) .* misFlag(:, poissFlag) ));

logLikeXbi = sum(sum( (log(muBi).*xBi + log((1-muBi)).*(1-xBi) ).*misFlag(:, binFlag) ));

logLikeX = logLikeXcont + logLikeXpo + logLikeXbi;


diffY = (y - eta);
logLikeY = sum( - log(sdy) - ( ((diffY).^2)./(2*(sdy.^2)) ) ); 

result = logLikeX + logLikeY;


function [remixedThetaStar] = remix(x, y, J, thetaStar)
 
global mu0 Sigma0 mu00 Sigma00 aBi bBi aPo bPo aSigmaY bSigmaY;
global abNuA abNuB sigBeta nCont nPo nBi nCov;
global contFlag poissFlag binFlag;
global eps;
global misFlagTrain;

u1 = sort(unique(J));
for i = 1:length(u1)
    X = x(J == i, :);
    Xcont = X(:, contFlag, :);
    Xpo = X(:, poissFlag, :);
    Xbi = X(:, binFlag, :);
    Xflagcont = misFlagTrain(J == i, contFlag, :);
    Xflagpo = misFlagTrain(J == i, poissFlag, :);
    Xflagbi = misFlagTrain(J == i, binFlag, :);
    
    Y = y(J == i, :);
    nY = length(J(J==i));
    nX = nY;
    
    mu = thetaStar(i).mu;
    sd = thetaStar(i).sd;
    muPo = thetaStar(i).muPo;
    muBi = thetaStar(i).muBi;    
    sdy = thetaStar(i).sdy;
    nuA = thetaStar(i).nuA;
    nuB = thetaStar(i).nuB;
    sigComp = thetaStar(i).sigComp;
    beta=thetaStar(i).beta;

    for j = 1:nCont
        mu(j) = getMu0(Xcont(Xflagcont(:, j)==1, j), mu(j), mu0(j), Sigma0(j), sd(j));   
        sd(j) = sqrt(getSig0(X(Xflagcont(:, j)==1, j), sd(j).^2, mu00(j), Sigma00(j), mu(j)));
    end

    
    for j = 1:nPo
        muPo(j) = gamrnd(aPo(j)+sum(Xpo(Xflagpo(:, j)==1, j)), 1./(1./bPo(j)+ nX) );
    end
    
    for j = 1:nBi
        sumX = sum(Xbi(Xflagbi(:, j)==1, j));
        muBi(j) = betarnd(aBi(j)+sumX, bBi(j)+ (nX - sumX) );
    end

        
    sdy = sqrt(getSig0(Y - [ones(nY, 1), X]*beta, sdy.^2, aSigmaY, bSigmaY, 0));
    
%     relatedBeta = beta(1, :);
%     nuA = sqrt(getSigBeta(relatedBeta, nuA^2, abNuA(1), abNuA(2)));
%     relatedBeta = reshape(beta(2:end, :), 1, nCov);
%     for repHype = 1:50
%     nuB = sqrt(getSigBeta(relatedBeta, nuB^2, abNuB(1), abNuB(2)));
%     end
%     sigComp = sigBeta; 
%     sigComp(1, :) = nuA*sigComp(1, :);
%     sigComp(2:end, :) = nuB*sigComp(2:end, :);


    [sampBeta, beta] = shortCut(sigComp, beta, X, Y, nY, sdy);   
    
    remixedThetaStar(i).mu = mu;
    remixedThetaStar(i).beta = beta;
    remixedThetaStar(i).sd = sd;
    remixedThetaStar(i).muPo = muPo;
    remixedThetaStar(i).muBi = muBi;    
    remixedThetaStar(i).sdy = sdy;
    remixedThetaStar(i).nuA = nuA;
    remixedThetaStar(i).nuB = nuB;
    remixedThetaStar(i).sigComp = sigComp;
    
    remixedThetaStar(i).sampBeta = sampBeta;

end



function [predY, postX, q] = getPredProb(x, thetaStar, nj, alpha)
global contFlag poissFlag binFlag nCont nPo nBi;
global misFlagTest;
global n nTest;

xCont = x(:, contFlag, :);
xPo = x(:, poissFlag, :);
xBi = x(:, binFlag, :);

iStar = length(nj);
q = zeros(nTest, iStar);

% This part is sample from G_0, which is used to get the predictive
% probability
% sd0 = repmat(sqrt(exp(normrnd(mu00, Sigma00))), nTest, 1);
% mu0 = repmat(normrnd(mu0, Sigma0), nTest, 1);
% nuA0 = sqrt(exp(normrnd(abNuA(1), abNuA(2))));
% nuB0 = sqrt(exp(normrnd(abNuB(1), abNuB(2))));
% sigComp0 = sigBeta;
% sigComp0(1, :) = nuA0*sigBeta(1, :);
% sigComp0(2:end, :) = nuB0*sigBeta(2:end, :);
% beta0 = normrnd(0, sigComp0);
% thetaStar(iStar+1).beta = beta0;

% This part uses the unique component to get the predictive probability
for k = 1:iStar
    mu = repmat(thetaStar(k).mu, nTest, 1);
    sd = repmat(thetaStar(k).sd, nTest, 1);
    sd = sd + .001;

    muBi = repmat(thetaStar(k).muBi, nTest, 1);
    muPo = repmat(thetaStar(k).muPo, nTest, 1);
    
    logLikeXcont = 0; logLikeXpo = 0; logLikeXbi = 0;

    if nCont > 0
    diffX = (xCont - mu);
    logLikeXcont = sum( - log(sd) - ( ((diffX).^2)./(2*(sd.^2)) ) .* misFlagTest(:, contFlag) , 2); 
    end
    
    if nPo > 0
    logLikeXpo = sum( ( log(muPo).*xPo - muPo ) .* misFlagTest(:, poissFlag) , 2);
    end
    
    if nBi >0
    logLikeXbi = sum( ( log(muBi).*xBi + log((1-muBi)).*(1-xBi) ) .* misFlagTest(:, binFlag) , 2);
    end
    
    q(:, k) = logLikeXcont + logLikeXpo + logLikeXbi;
    q(:, k) = q(:, k) + log(nj(k)) - log(n+alpha);    
end

% diffX = (x - mu0);
% q(:, iStar+1) = sum( - log(sd0) - ( ((diffX).^2)./(2*(sd0.^2)) ), 2);     
% q(:, iStar+1) = q(:, iStar+1) + log(alpha) - log(n+alpha);

% This calculates P(x) for the current iteration
m = max(q');
postX(:, 1) = exp(m'+log(sum(exp(q - repmat(m', 1, iStar) ), 2)));


% This part gets P(y) for the current iteration
qRel = q - repmat(m', 1, iStar);
q = exp(qRel);
q = q./repmat(sum(q, 2), 1, iStar);

for k = 1:iStar
    etaT = [ones(nTest, 1), x]*thetaStar(k).beta;
    predY(:, k) = etaT.*q(:, k);
end
    
predY = sum(predY, 2);



function [sampBeta, beta] = shortCut(sigComp, beta, X, Y, nY, sdy)
global epsVec;

e = ( 1 ./ sqrt( 1 ./ (sigComp .^ 2) + nY/4) );

M = 2;
L = 2;

E = getE([ones(nY, 1), X], Y, sdy, beta, sigComp);
g = getG([ones(nY, 1), X], Y, sdy, beta, sigComp);

sampBeta = [];
for contEps = 1:length(epsVec)
    e = epsVec(contEps)*e;
    direction = 1; bounce = 0; newF = []; newB = [];
    m = 1;
    while m <= M & direction == 1
        for l = 1:L            
            [beta, E, g, accept(l)] = getBeta([ones(nY, 1), X], Y, sdy, beta, sigComp, E, g, e);
        end

        lenF = length(newF); 
        newF(lenF+1).beta = beta;
        newF(lenF+1).E = E;
        newF(lenF+1).g = g;
        newF(lenF+1).weight = 1;   
        lenF = length(newF);
        m = m+1;
        if sum(accept)==0 
            direction = -1;
            ind = lenF;
            inc = -1;
            ind = ind+inc;                    
            while m<=M & ind > 0
                newF(ind).weight = 2;
                beta = newF(ind).beta;
                E = newF(ind).E;
                g = newF(ind).g;                        
                ind = ind+inc;
                m = m+1;
            end
        end
    end

    
    while m <= M & direction == -1
        for l = 1:L
            [beta, E, g, accept(l)] = getBeta([ones(nY, 1), X], Y, sdy, beta, sigComp, E, g, e);
        end

        lenB = length(newB); 
        newB(lenB+1).beta = beta;
        newB(lenB+1).E = E;
        newB(lenB+1).g = g;
        newB(lenB+1).weight = 1;   
        lenB = length(newB);
        m = m+1;
        if sum(accept)==0 
            bounce = 1;
            break;
        end
    end
    
   
    if m<=M & bounce == 1
        lenB = length(newB);        
        lenF = length(newF);  
        lenS = length(sampBeta);
        for temp1 = 1:lenB
            sampBeta(lenS+temp1).beta = newB(temp1).beta; 
            sampBeta(lenS+temp1).weight = newB(temp1).weight; 
        end
        lenS = length(sampBeta);        
        for temp1 = 1:lenF
            sampBeta(lenS+temp1).beta = newF(temp1).beta; 
            sampBeta(lenS+temp1).weight = newF(temp1).weight; 
        end
        lenS = length(sampBeta); 
        ind = 1;
        while m<=M
            if ind == 1
                inc = 1;
            elseif ind == lenS; 
                inc = -1;
            end
            sampBeta(ind).weight = sampBeta(ind).weight+1;
            beta = sampBeta(ind).beta;
            m = m+1;
            ind = ind+inc;
        end        
    else
        lenB = length(newB);        
        lenF = length(newF);  
        lenS = length(sampBeta);  
        for temp1 = 1:lenB
            sampBeta(lenS+temp1).beta = newB(temp1).beta; 
            sampBeta(lenS+temp1).weight = newB(temp1).weight; 
        end
        lenS = length(sampBeta);        
        for temp1 = 1:lenF
            sampBeta(lenS+temp1).beta = newF(temp1).beta; 
            sampBeta(lenS+temp1).weight = newF(temp1).weight; 
        end                
    end
            
end

lenS = length(sampBeta);



function [postX, q] = getPostX(x, thetaStar, nj, alpha)
global contFlag poissFlag binFlag nCont nPo nBi;
global misFlagTest;
global n nTest;

xCont = x(:, contFlag, :);
xPo = x(:, poissFlag, :);
xBi = x(:, binFlag, :);

iStar = length(nj);
q = zeros(nTest, iStar);

% This part is sample from G_0, which is used to get the predictive
% probability
% sd0 = repmat(sqrt(exp(normrnd(mu00, Sigma00))), nTest, 1);
% mu0 = repmat(normrnd(mu0, Sigma0), nTest, 1);
% nuA0 = sqrt(exp(normrnd(abNuA(1), abNuA(2))));
% nuB0 = sqrt(exp(normrnd(abNuB(1), abNuB(2))));
% sigComp0 = sigBeta;
% sigComp0(1, :) = nuA0*sigBeta(1, :);
% sigComp0(2:end, :) = nuB0*sigBeta(2:end, :);
% beta0 = normrnd(0, sigComp0);
% thetaStar(iStar+1).beta = beta0;

% This part uses the unique component to get the predictive probability
for k = 1:iStar
    mu = repmat(thetaStar(k).mu, nTest, 1);
    sd = repmat(thetaStar(k).sd, nTest, 1);
    sd = sd + .001;

    muBi = repmat(thetaStar(k).muBi, nTest, 1);
    muPo = repmat(thetaStar(k).muPo, nTest, 1);
    
    logLikeXcont = 0; logLikeXpo = 0; logLikeXbi = 0;

    if nCont > 0
    diffX = (xCont - mu);
    logLikeXcont = sum( - log(sd) - ( ((diffX).^2)./(2*(sd.^2)) ) .* misFlagTest(:, contFlag) , 2); 
    end
    
    if nPo > 0
    logLikeXpo = sum( ( log(muPo).*xPo - muPo ) .* misFlagTest(:, poissFlag) , 2);
    end
    
    if nBi >0
    logLikeXbi = sum( ( log(muBi).*xBi + log((1-muBi)).*(1-xBi) ) .* misFlagTest(:, binFlag) , 2);
    end
    
    q(:, k) = logLikeXcont + logLikeXpo + logLikeXbi;
    q(:, k) = q(:, k) + log(nj(k)) - log(n+alpha);    
end

% diffX = (x - mu0);
% q(:, iStar+1) = sum( - log(sd0) - ( ((diffX).^2)./(2*(sd0.^2)) ), 2);     
% q(:, iStar+1) = q(:, iStar+1) + log(alpha) - log(n+alpha);

% This calculates P(x) for the current iteration
m = max(q');
postX(:, 1) = exp(m'+log(sum(exp(q - repmat(m', 1, iStar) ), 2)));


% This part gets P(y) for the current iteration
qRel = q - repmat(m', 1, iStar);
q = exp(qRel);
q = q./repmat(sum(q, 2), 1, iStar);

function postY = getPostY(x, beta)
postY = [ones(size(x, 1), 1), x]*beta;




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

logLike = (alpha^iStar)*(exp(gammaln(alpha)-gammaln(alpha+n))); 

diffSig = x - a0;
logPrior = sum(-.5*log(2*pi)-log(b0)-( (( diffSig ).^2)./(2*(b0.^2)) ));

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


function [updateB, updateE, updateG, accept] = getBeta(x, y, sdy, b, sigma0, E, g, e)
global leapFrog;

[dim1, dim2] = size(b);

oldB = b;
oldE = E;
oldG = g;
if E == 0;
    E = getE(x, y, sdy, oldB, sigma0);
    g = getG(x, y, sdy, oldB, sigma0);
end

    
    p = randn(size(b));
    H = sum(sum( .5*(p.^2) ))+ E;
    newB = b; 
    newG = g; 

    for leap=1:leapFrog
        p = p - e.*newG/2;
        newB = newB + e.*p;
        tempG = getG(x, y, sdy, newB, sigma0);
        if isfinite(tempG)
            newG = tempG;
        end

        p = p - e.*newG/2;
    end

    newE = getE(x, y, sdy, newB, sigma0);
    newH = sum(sum( .5*(p.^2) )) + newE;

    acceptProb = min(1, exp(H - newH));

    if (rand < acceptProb)
        updateB = newB;
        updateE = newE;
        updateG = newG;
        accept = 1;
    else
        updateB = oldB;
        updateE = oldE;
        updateG = oldG;
        accept = 0;
    end
   

% This calculates the energy function for the posterior distribution. 

function E = getE(x, y, sdy, beta, sigma)
global nResp;

[n, d] = size(y);

eta = x*beta;
sdy = repmat(sdy, n, 1);

diffY = (y - eta);
logLike = sum(sum( - log(sdy) - ( ((diffY).^2)./(2*(sdy.^2)) ) )); 
logPrior =  sum(sum( ( -0.5*(beta).^2 ) ./ (sigma.^2) ));

E = -(logLike + logPrior);



% This part calculates the derivatives for all parameters. 
% Note that the calculation is completely vectorized and quite fast. 
% Moreover, the code is written in a way that avoids overflow.

function g = getG(x, y, sdy, beta, sigma)
global nResp;

[n, d] = size(y);

eta = x*beta;
sdy = repmat(sdy, n, 1);

dLogLike =  x'*( (y - eta) ./ (sdy.^2) );

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



function createfigure(XData1, YData1, XData2, YData2, X1, Y1, X2, Y2)
%CREATEFIGURE(XDATA1,YDATA1,XDATA2,YDATA2,X1,Y1,X2,Y2)
%  XDATA1:  line xdata
%  YDATA1:  line ydata
%  XDATA2:  line xdata
%  YDATA2:  line ydata
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 27-Jul-2007 09:22:05

% Create figure
figure1 = figure('XVisual',...
    '0x22 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
axes1 = axes('Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
% xlim([-2.637 2.5]);
hold('all');

% Create line
line(XData1,YData1,'Parent',axes1,'MarkerSize',15,'Marker','.',...
    'LineStyle','none',...
    'Color',[1 0 0]);

% Create line
line(XData2,YData2,'Parent',axes1,'MarkerSize',15,'Marker','.',...
    'LineStyle','none',...
    'Color',[0 1 1]);

% Create xlabel
xlabel('Dem 1980','FontSize',14);

% Create ylabel
ylabel('Dem 1984','FontSize',14);

% Create line
line(XData1,YData1,'Parent',axes1,'MarkerSize',15,'Marker','.',...
    'LineStyle','none');

% Create line
line(XData2,YData2,'Parent',axes1,'MarkerSize',15,'Marker','x',...
    'LineStyle','none');

% Create plot
plot(X1,Y1,'Parent',axes1,'Color',[0 0 0]);

% Create plot
plot(X2,Y2,'Parent',axes1,'Color',[0 0 0]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.768 0.8064 0.125 0.09359],'FontSize',12);


function [xTrain, yTrain, xTest, yTest, zTrain, zTest] = createData(n, nTrain, state)

rand('state', state);
randn('state', state);

    
x1 = normrnd(-1, 1, n/2, 1);
x2 = normrnd(1, 1, n/2, 1);

a1 = -1; a2 = 1;
% b = [normrnd(1, 1, nCov, 1)];
b = 0.5;

y1 = [ones(n/2, 1), x1]*[a1; b] + normrnd(0, .25, n/2, 1);
y2 = [ones(n/2, 1), x2]*[a2; b] + normrnd(0, .25, n/2, 1);

x = [x1; x2];

y = [y1; y2];

z = [-1*ones(n/2, 1); ones(n/2, 1)];

trainInd = randsample(n, nTrain);
testInd  = setdiff([1:n], trainInd);

data = [y, x];

train = data(trainInd, :);
test  = data(testInd, :);

xTrain = train(:, 2:end);
yTrain = train(:, 1);
xTest  = test(:, 2:end);
yTest  = test(:, 1);

zTrain = z(trainInd, :);
zTest = z(testInd, :);

