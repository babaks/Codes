% This file imports the 50 datasets for simulation 1 and uses the DPM model
% to analyze each dataset. The main output here is "nu", which
% as described in the paper is a measure of significance (i.e., the higher
% nu the more significant is the variable). The variables with nu close to 
% zero or negative are basically dispensable. 

function nu = dpmSim2(numberOfIterations, burnIn)
global contFlag poissFlag binFlag 
n = 30; % Sample size
% These are the flags to identify continuous, categorical and count variables. 
contFlag = []; 
binFlag = [1:3];  
poissFlag = [];
    
impData = importdata('Simulation2.txt');    
allData = impData.data;

for i = 1:50
    i
    data = allData(allData(:, 1)==i, 2:end);
    thisSamp = DPM(data, numberOfIterations, burnIn);

    % devRepAvg is the expectected deviance for each replicated dataset,
    % devObsAvg is the expected deviance based for the observed data and
    % pDobs is the penalty added to the observed deviance. 
    nu(i, :) = mean(thisSamp.devRepAvg, 1) - (thisSamp.devObsAvg + thisSamp.pDobs);
end


function thisSamp = DPM(data, numberOfIterations, burnIn)
global nTest n nBi nLeaf nCov;
global contFlag poissFlag binFlag missFlagTest;
global mu0 Sigma0 mu00 Sigma00 aBi aPo bPo a0 b0;
global aSigma0 bSigma0 muMu sigMu;
global aSigma00 bSigma00 muSig sigSig;
global abNuA abNuB sigBeta;
global leapFrog epsVec;
global dummyInd;

x = data(:, 2:end);
y = data(:, 1);
yD= dummyvar(y);

n = size(x, 1);

nCov = size(x, 2);
source(1).ind = [1:nCov];
nLeaf = size(yD, 2);

nCont = length(contFlag);
nBi = length(binFlag);
nPo = length(poissFlag);

X.Cont = x(:, contFlag);
X.ContMiss = ones(size(X.Cont));
X.ContMiss(find(isnan(X.Cont))) = 0; 
X.Cont(find(isnan(X.Cont))) = 0;

X.Cont = X.Cont - repmat(nanmean(X.Cont), n, 1);
X.Cont = X.Cont ./ repmat(nanstd(X.Cont), n, 1);

X.Po = x(:, poissFlag);
X.PoMiss = ones(size(X.Po));
X.PoMiss(find(isnan(X.Po))) = 0; 
X.Po(find(isnan(X.Po))) = 0;

X.Bi = x(:, binFlag);
X.BiMiss = ones(size(X.Bi));

for i = 1:nBi
    ind = find(~isnan(X.Bi(:, i)));
    X.BiDummy(i).a(ind, :) = dummyvar(X.Bi(ind, i));
end

m = max(X.Bi);
dummyInd = cumsum([1, m]);

X.BiMiss(find(isnan(X.Bi))) = 0; 
X.Bi(find(isnan(X.Bi))) = 0; 

biDummy =  cat(2, X.BiDummy.a);
for i = 1:n
   
Z(i, 1) = struct('cont', X.Cont(i, :), 'bi', X.Bi(i, :), 'biDummy', biDummy(i, :), 'po', X.Po(i, :), ...
    'contMiss', X.ContMiss(i, :), 'biMiss', X.BiMiss(i, :), 'poMiss', X.PoMiss(i, :));
end

XT.Cont = [];
XT.Bi = [];
XT.Po = [];

xTest = zeros(n, nCov);
xTest(:, contFlag) = X.Cont;
xTest(:, binFlag) = X.Bi;
xTest(:, poissFlag) = X.Po;
xTest = repmat(xTest, nCov, 1);
nTest = size(xTest, 1);

epsVec = [.01, 0.5, 1];  

leapFrog = 50;
nSkip = 3;


a0 = -2; b0 = 3; 
mu0 = zeros(1, nCont);
Sigma0 = 5*ones(1, nCont);
mu00 = zeros(1, nCont);
Sigma00 = 2*ones(1, nCont);

aSigma00 = -1; bSigma00 = 2;

for i = 1:nBi
    aBi(i).a = 1*ones(1, max(X.Bi(:, i)));
end

aPo = .5*ones(1, nPo); bPo = 10*ones(1, nPo);

abSigmaInt = [-1, 2];
abSigma = [-1, 2];
abEta = [-2, 3];
abNuA  = [0, .1];
abNuB  = [0, .1];

intSigma = 1;
ardSigma = ones(nCov, nLeaf);
eta = ones(1, length(source));

nuA = 1;
nuB = 1;
sigBeta = ones(nCov+1, nLeaf);
sigBeta(1, :) = intSigma;
sigBeta(2:end, :) = ardSigma;

sigComp = sigBeta;
sigComp(1, :) = 5*nuA*sigBeta(1, :);
sigComp(2:end, :) = 5*nuB*sigBeta(2:end, :);

alpha(1) = .001;
thetaStar.mu = nanmean(x(:, contFlag));
thetaStar.sd = nanstd(x(:, contFlag));
thetaStar.muPo = 1*ones(1, nPo);
for i = 1:nBi
thetaStar.muBi(i).a = aBi(i).a / sum(aBi(i).a);
end
thetaStar.nuA = nuA;
thetaStar.nuB = nuB;
thetaStar.sigComp = sigComp;
thetaStar.beta = ones(1+nCov, nLeaf);

newJ = repmat(1, n, 1);
newNj = n;

countP = 0; 
for iter = 1:numberOfIterations
iter


[thetaStar, newJ, newNj] = main_mvn(Z, yD, thetaStar, newJ, newNj, alpha(iter));
newNj'
nComp = length(newNj);

thetaStar = remix_mvn(Z, yD, newJ, thetaStar);

muT1 = cat(1, thetaStar.mu);
sdT1 = cat(1, thetaStar.sd);

for j = 1:nCont
    muT2 = muT1(newJ, j);
    sdT2 = sdT1(newJ, j);
    
    XT.Cont(:, j) = normrnd(muT2, sdT2);
    
end

for i = 1:nComp
for j = 1:nBi
    muT1 = thetaStar(i).muBi(j).a;
    temp = mnrnd(1, muT1, newNj(i));    
    [m, XT.Bi(newJ==i, j)] = max(temp, [], 2);
end
end

muT1 = cat(1, thetaStar.muPo);
for j = 1:nPo
    muT2 = muT1(newJ, j);
    
    XT.Po(:, j) = poissrnd(muT2);
    
end



alpha(iter+1) = pickAlpha(alpha(iter), length(newNj), a0, b0);

nComp = length(newNj);

if iter >= burnIn
    if rem(iter, nSkip) == 0
    countP = countP+1;
    
    
    samp(countP).thetaStar = thetaStar;
    samp(countP).newNj = newNj;
    samp(countP).J = newJ;    
    samp(countP).alpha = alpha(iter+1);
    samp(countP).sigBeta = sigBeta;
    samp(countP).XT = XT;

    end
end



if iter >= 5

    
    interceptBeta = [];
    for k = 1:nComp
        interceptBeta = cat(2, interceptBeta, thetaStar(k).beta(1, :)/thetaStar(k).nuA);
    end
    
    for repHype = 1:3
        intSigma = sqrt(getSigBeta(interceptBeta, intSigma^2, abSigmaInt(1), abSigmaInt(2)));
    end
    sigBeta(1, :) = intSigma;


    newBeta = [];
    for k = 1:nComp
        tempB = thetaStar(k).beta(2:end, :);
        newBeta = [newBeta, tempB];
    end
       
    tempSigma = ardSigma;        
    for s = 1:length(source)
        relatedBeta = newBeta(source(s).ind, :);
        reshapedBeta = reshape(relatedBeta, 1, (length(source(s).ind))*nLeaf*nComp);  % Please note, eta does not control the last covariate          
        for repHype = 1:10
            eta(s) = sqrt(getSigBeta(reshapedBeta, eta(s)^2, abEta(1), abEta(2)));        
        end
    end
    
    sigBeta(2:end, :) = eta;


end
 

end

thisSamp.J = cat(2, samp.J);

nSamp = length(samp);
predProbObs = [];
missFlagTest = [];
missFlagTest(:, cat(2, contFlag, binFlag, poissFlag)) = cat(2, X.ContMiss, X.BiMiss, X.PoMiss);
xTrain = [];
xTrain(:, cat(2, contFlag, binFlag, poissFlag)) = cat(2, X.Cont, X.Bi, X.Po);



sumProb = 0;
devObs = [];
allWeights = [];
nSampP = 0;
count = 0;
for j = 1:1:nSamp
    
    count = count+1;

    xTrain2 = xTrain;
    xTrain2(:, binFlag) = xTrain2(:, binFlag)-1;
        
    thetaStar = samp(j).thetaStar;
    newNj = samp(j).newNj;
    alpha = samp(j).alpha;
    sigBeta = samp(j).sigBeta;

    [pX, q] = getPostX(xTrain, thetaStar, newNj); 
    pY = [];
    thisP = [];
    for temp1 = 1:length(thetaStar)
        thisBeta = thetaStar(temp1).beta;
        thisP = getPostY(xTrain2, thisBeta).*repmat(q(:, temp1), 1, nLeaf);
        pY = cat(3, pY, thisP);
    end

    pY = sum(pY, 3);
    
    devObs = cat(1, devObs, (-2*sum((sum(yD.*log(pY), 2)), 1)) );
    
    sumProb = (pY+sumProb);    
    nSampP = nSampP + 1;
   
end

predProbObs = sumProb(:, 2)/nSampP;
devObsAvg = mean(devObs);
pDobs = 0.5 * (1/(length(devObs)-1))* sum( ((devObs - devObsAvg).^2));

thisSamp.predProbObs = predProbObs;
thisSamp.devObsAvg = devObsAvg;
thisSamp.pDobs = pDobs;


nSamp = length(samp);
predProb = [];
devRepAvg = [];
pDrep = [];

missFlagTest = [];
missFlagTest(:, cat(2, contFlag, binFlag, poissFlag)) = cat(2, X.ContMiss, X.BiMiss, X.PoMiss);
missFlagTest = repmat(missFlagTest, nCov, 1);
for i = 1:nSkip:nSamp
    thisTest(:, cat(2, contFlag, binFlag, poissFlag)) = cat(2, samp(i).XT.Cont, samp(i).XT.Bi, samp(i).XT.Po);
    for j = 1:nCov
    xTest(((j-1)*n+1):(j*n), j) = thisTest(:, j);
    end


    xTest2 = xTest;
    xTest2(:, binFlag) = xTest2(:, binFlag)-1;

    
    
    sumProb = 0;
    devRep = [];
    nSampP = 0;
    count = 0;

for j = 1:1:nSamp;
    count = count+1;
    thetaStar = samp(j).thetaStar;
    newNj = samp(j).newNj;
    alpha = samp(j).alpha;
    sigBeta = samp(j).sigBeta;

    [pX, q] = getPostX(xTest, thetaStar, newNj); 
    pY = [];
    for temp1 = 1:length(thetaStar)
        thisBeta = thetaStar(temp1).beta;
        thisP = getPostY(xTest2, thisBeta).*repmat(q(:, temp1), 1, nLeaf);
            
        pY = cat(3, pY, thisP);
            
    end

    pY = sum(pY, 3);
    
    thisDevRep = [];
    for k = 1:nCov
        thisDevRep(:, k) = -2*sum( (sum(yD.* log( pY(((k-1)*n+1):(k*n), :) ), 2)), 1);
    end

    devRep = cat(1, devRep, thisDevRep);

    
    sumProb = (pY+sumProb);    
    nSampP = nSampP + 1;
end

thisAvg = mean(devRep, 1);
devRepAvg = cat(1, devRepAvg, thisAvg);
pDrep = cat(1, pDrep, 0.5 * (1/(size(devRep, 1)-1))* sum( ((devRep - repmat(thisAvg, size(devRep, 1), 1)).^2) ));


predProb = cat(2, predProb, sumProb(:, 2)/nSampP);
end
    
thisSamp.predProb = predProb;
thisSamp.devRepAvg = devRepAvg;
thisSamp.pDrep = pDrep;




function [thetaStar, updatedJ, updatedNj] = main_mvn(x, y, thetaStar, J, nj, alpha);

global mu0 Sigma0 mu00 Sigma00 aBi aPo bPo;
global n nBi sigBeta;

M=2;

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
            for ii = 1:nBi
            tempGam = gamrnd(aBi(ii).a, 1);
            muBi(ii).a = tempGam/sum(tempGam);
            end

            nuA = 1;
            nuB = 1;
            sigComp = sigBeta;
            sigComp(1, :) = nuA*sigBeta(1, :);
            sigComp(2:end, :) = nuB*sigBeta(2:end, :);
            beta = normrnd(0, sigComp);            
            thetaStar(kBar+1+m).mu = mu;
            thetaStar(kBar+1+m).sd = sd;
            thetaStar(kBar+1+m).muPo = muPo;
            thetaStar(kBar+1+m).muBi = muBi;
            
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
            for ii = 1:nBi
            tempGam = gamrnd(aBi(ii).a, 1);
            muBi(ii).a = tempGam/sum(tempGam);
            end
            nuA = 1;
            nuB = 1;
            sigComp = sigBeta;
            sigComp(1, :) = nuA*sigBeta(1, :);
            sigComp(2:end, :) = nuB*sigBeta(2:end, :);
            beta = normrnd(0, sigComp);            
            thetaStar(kBar+m).mu = mu;
            thetaStar(kBar+m).sd = sd;
            thetaStar(kBar+m).muPo = muPo;
            thetaStar(kBar+m).muBi = muBi;            
            
            thetaStar(kBar+m).nuA = nuA;
            thetaStar(kBar+m).nuB = nuB;
            thetaStar(kBar+m).sigComp = sigComp;
            thetaStar(kBar+m).beta = beta;
         end
    end

    q1 = zeros(kBar, 1);
        
    for k = 1:kBar
        mu = thetaStar(k).mu;
        sd = thetaStar(k).sd;        
        muPo = thetaStar(k).muPo;
        muBi = thetaStar(k).muBi;
        
        X = cat(2, x(i).cont, x(i).bi-1, x(i).po);
        
        beta = thetaStar(k).beta;
        eta = [1, X]*beta;
        sigComp = thetaStar(k).sigComp;
        q1(k) = mvnLogLike(x(i), y(i, :), mu, sd, muPo, muBi, eta);
    end
    q1 = q1 + log(nj) - log(n-1+alpha);
    
    q2 = zeros(M, 1);
    for k = 1:M
        mu = thetaStar(kBar+k).mu;
        sd = thetaStar(kBar+k).sd;
        muPo = thetaStar(kBar+k).muPo;
        muBi = thetaStar(kBar+k).muBi;
        
        X = cat(2, x(i).cont, x(i).bi-1, x(i).po);
        
        beta = thetaStar(kBar+k).beta;
        eta = [1, X]*beta;
        sigComp = thetaStar(kBar+k).sigComp;
        q2(k) = mvnLogLike(x(i), y(i, :), mu, sd, muPo, muBi, eta);
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



function result = mvnLogLike(x, y, mu, sd, muPo, muBi, eta)
global nLeaf nBi dummyInd;

nX = size(x, 1);

mu = repmat(mu, nX, 1);
sd = repmat(sd, nX, 1);
muPo = repmat(muPo, nX, 1);
muBi = repmat(muBi, nX, 1);

sd = sd + .001;
sd(sd==Inf)=max(max(sd(sd~=Inf)));
diffX = (x.cont - mu);
logLikeXcont = sum(sum( ( - log(sd) - ( ((diffX).^2)./(2*(sd.^2)) ) ) .* x.contMiss)); 


logLikeXpo = sum(sum( ( log(muPo).*x.po - muPo ) .* x.poMiss ));

for i = 1:nBi
    logLikeXbi(i) = log(mnpdf(x.biDummy(dummyInd(i):dummyInd(i+1)-1), muBi(i).a));
end
logLikeXbi = sum(logLikeXbi);

m = max(eta');
logLikeY = sum( (sum((y.*eta), 2) - (m' + log(sum(exp(eta - repmat(m', 1, nLeaf) ), 2)) ) ) );

result = sum([logLikeXcont, logLikeXbi, logLikeXpo, logLikeY], 2);



function remixedThetaStar = remix_mvn(x, y, J, thetaStar)
 
global mu0 Sigma0 mu00 Sigma00 aPo bPo aBi;
global nCont nPo nBi;
global dummyInd;

u1 = sort(unique(J));
for i = 1:length(u1)
    nX = length(find(J==i));
    X = x(J == i);
    
    mu = thetaStar(i).mu;
    sd = thetaStar(i).sd;
    muPo = thetaStar(i).muPo;
    muBi = thetaStar(i).muBi;
    
    Y = y(J == i, :);
    nY = length(J(J==i));

    nuA = thetaStar(i).nuA;
    nuB = thetaStar(i).nuB;
    sigComp = thetaStar(i).sigComp;
    beta=thetaStar(i).beta;

    xCont = cat(1, X.cont);
    contMiss = cat(1, X.contMiss);
    for j = 1:nCont
        mu(j) = getMu0(xCont(contMiss(:, j)==1, j), mu(j), mu0(j), Sigma0(j), sd(j));   
        sd(j) = sqrt(getSig0(xCont(contMiss(:, j)==1, j), sd(j).^2, mu00(j), Sigma00(j), mu(j)));
    end

    xPo = cat(1, X.po);
    poMiss = cat(1, X.poMiss);
    for j = 1:nPo
        muPo(j) = gamrnd(aPo(j)+sum(xPo(poMiss(:, j)==1, j)), 1./(1./bPo(j)+ nX) );
    end

    
    biDummy = cat(1, X.biDummy);
    for j = 1:nBi
         thisX = (biDummy(:, dummyInd(j):dummyInd(j+1)-1));
         sumX = sum(thisX, 1);
         temp = (aBi(j).a+sumX);
         tempGam = gamrnd(temp, 1);
         muBi(j).a = tempGam/sum(tempGam);
    end
    

    thisX = cat(2, cat(1, X.cont), cat(1, X.bi)-1, cat(1, X.po));
    [sampBeta, beta] = shortCut(sigComp, beta, thisX, Y, nY);
    
    remixedThetaStar(i).mu = mu;
    remixedThetaStar(i).sd = sd;
    remixedThetaStar(i).muPo = muPo;
    remixedThetaStar(i).muBi = muBi;    
    remixedThetaStar(i).beta = beta;
    remixedThetaStar(i).nuA = nuA;
    remixedThetaStar(i).nuB = nuB;
    remixedThetaStar(i).sigComp = sigComp;
    
    remixedThetaStar(i).sampBeta = sampBeta;
    
end


function [sampBeta, beta] = shortCut(sigComp, beta, X, Y, nY)
global epsVec;

e = ( 1 ./ sqrt( 1 ./ (sigComp .^ 2) + nY/4) );

M = 2;
L = 2;

E = getE([ones(nY, 1), X], Y, beta, sigComp);
g = getG([ones(nY, 1), X], Y, beta, sigComp);

sampBeta = [];
for contEps = 1:length(epsVec)
    e = epsVec(contEps)*e;
    direction = 1; bounce = 0; newF = []; newB = [];
    m = 1;
    while m <= M && direction == 1
        for l = 1:L
            [beta, E, g, accept(l)] = getBeta([ones(nY, 1), X], Y, beta, sigComp, E, g, e);
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
            [beta, E, g, accept(l)] = getBeta([ones(nY, 1), X], Y, beta, sigComp, E, g, e);
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


function [postX, q] = getPostX(z, thetaStar, nj)
global n missFlagTest;
global nCont nPo nBi;
global contFlag poissFlag binFlag;

nTest = size(z, 1);

x.cont = z(:, contFlag);
x.po = z(:, poissFlag);
x.bi = z(:, binFlag);

x.contMiss = missFlagTest(:, contFlag);
x.poMiss = missFlagTest(:, poissFlag);
x.biMiss = missFlagTest(:, binFlag);

iStar = length(nj);
q = zeros(nTest, iStar);

for k = 1:iStar

    mu = repmat(thetaStar(k).mu, nTest, 1);
    sd = repmat(thetaStar(k).sd, nTest, 1);
    muPo = repmat(thetaStar(k).muPo, nTest, 1);
    muBi = thetaStar(k).muBi;
    
sd = sd + .001;
sd(sd==Inf)=max(max(sd(sd~=Inf)));
diffX = (x.cont - mu);
logLikeXcont = (( ( - log(sd) - ( ((diffX).^2)./(2*(sd.^2)) ) ) .* x.contMiss)); 


logLikeXpo = (( ( log(muPo).*x.po - muPo ) .* x.poMiss ));

logLikeXbi = zeros(nTest, nBi);
for i = 1:nBi
    logLikeXbi(x.biMiss(:, i)==1, i) = (log(sum(repmat(muBi(i).a, sum(x.biMiss(:, i)==1), 1) .* dummyvar(x.bi(x.biMiss(:, i)==1, i)), 2)));
end

    q(:, k) = sum(cat(2, logLikeXcont, logLikeXpo, logLikeXbi), 2);
    q(:, k) = q(:, k) + log(nj(k)) - log(n);    
end


m = max(q');
postX(:, 1) = exp(m'+log(sum(exp(q - repmat(m', 1, iStar) ), 2)));

qRel = q - repmat(m', 1, iStar);
q = exp(qRel);
q = q./repmat(sum(q, 2), 1, iStar);


function postY = getPostY(x, beta)
global nLeaf;

etaT = [ones(size(x, 1), 1), x]*beta;
m = max(etaT');        
postY = exp(etaT - repmat(m' + log(sum(exp(etaT - repmat(m', 1, nLeaf) ), 2)), [1, nLeaf, 1] ) );


function newAlpha = pickAlpha(oldAlpha, iStar, a0, b0)
m = 20; w = 2;

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




function newMu = getMu0(y, oldMu, mu0, sigma0, sigma)
m = 30; w = 2;

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
m = 20; w = 2;

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



function [updateB, updateE, updateG, accept] = getBeta(d, r, b, sigma0, E, g, e)
global leapFrog;

[dim1, dim2] = size(b);

oldB = b;
oldE = E;
oldG = g;
if E == 0;
    E = getE(d, r, oldB, sigma0);
    g = getG(d, r, oldB, sigma0);
end
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
        updateE = newE;
        updateG = newG;
        accept = 1;
    else
        updateB = oldB;
        updateE = oldE;
        updateG = oldG;
        accept = 0;
    end

function E = getE(d, r, beta, sigma)
global nLeaf;

eta = d*beta;

m = max(eta');
logLike = sum( (sum((r.*eta), 2) - (m' + log(sum(exp(eta - repmat(m', 1, nLeaf) ), 2)) ) ) );

logPrior =  sum(sum( ( -0.5*(beta).^2 ) ./ (sigma.^2) ));

E = -(logLike + logPrior);


function g = getG(d, r, beta, sigma)
global nLeaf;
eta = d*beta;

m = max(eta');

dLogLike = (d'*r - d'*exp( eta - repmat( (m' + log(sum(exp(eta - repmat(m', 1, nLeaf) ), 2)) ) , 1, nLeaf) ) );
dLogPrior =  -(beta) ./ (sigma.^2);
 
g = -(dLogLike + dLogPrior);





function newSigma = getSigBeta(beta, oldSigma2, mu, tau)
m = 20; w = 2;

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


