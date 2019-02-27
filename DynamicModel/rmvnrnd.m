function [X, rho, nar, ngibbs] = rmvnrnd(mu,sigma,N,A,b,rhoThr,debug)
%RMVNRND Draw from the  truncated multivariate normal distribution.
%   X = rmvnrnd(MU,SIG,N,A,B) returns in N-by-P matrix X a
%   random sample drawn from the P-dimensional multivariate normal
%   distribution with mean MU and covariance SIG truncated to a
%   region bounded by the hyperplanes defined by the inequalities Ax<=B.
%
%   [X,RHO,NAR,NGIBBS]  = rmvnrnd(MU,SIG,N,A,B) returns the
%   acceptance rate RHO of the accept-reject portion of the algorithm
%   (see below), the number NAR of returned samples generated by
%   the accept-reject algorithm, and the number NGIBBS returned by
%   the Gibbs sampler portion of the algorithm.
%
%   rmvnrnd(MU,SIG,N,A,B,RHOTHR) sets the minimum acceptable
%   acceptance rate for the accept-reject portion of the algorithm
%   to RHOTHR. The default is the empirically identified value
%   2.9e-4.
%
%   ALGORITHM
%   The function begins by drawing samples from the untruncated MVN
%   distribution and rejecting those which fail to satisfy the
%   constraints. If, after a number of iterations with escalating
%   sample sizes, the acceptance rate is less than RHOTHR it
%   switches to a Gibbs sampler. 
%
% ACKNOWLEDGEMENT
%   This makes use of TruncatedGaussian by Bruno Luong (File ID:
%   #23832) to generate draws from the one dimensional truncated normal.
%
%   REFERENCES
%   Robert, C.P, "Simulation of truncated normal variables",
%   Statistics and Computing, pp. 121-125 (1995).




% Copyright 2011-2013 Tim J. Benham, School of Mathematics and Physics,
%                University of Queensland.

%if debug, disp('rmvnrnd'), end

%
% Constant parameters
%
defaultRhoThr = 2.9e-4;                   % min. acceptance rate to
                                          % accept accept-reject sampling.
%
% Process input arguments.
%
if nargin<6 || rhoThr<0,rhoThr = defaultRhoThr;end

if exist('A', 'var')
    A = A';
else
    A = [];
end
if exist('b', 'var')
    b = b';
else
    b = [];
end
if exist('mu', 'var')
    p = length(mu);                        % dimension -- mu can't
                                           % be defaulted!
    if size(mu,2) == 1
        mu = mu';
    end
else
    error('mean vector mu must be supplied')
end

m = size(A,2);                          % number of constraints
if m==0
    A = zeros(p,1);
    b = zeros(1);
end

% initialize return arguments
X = zeros(N,p);                         % returned random variates
nar = 0;                                % no. of X generated by a-r
ngibbs = 0;
rho = 1; 

if rhoThr<1
    % Approach 1: accept-reject method
    n = 0; % no. accepted
    maxSample = 1e6;
    trials = 0; passes = 0;
    s = N;
    while n<N && ( rho>rhoThr || s<maxSample)
        R  = mvnrnd(mu,sigma,s);
        R = R(sum(R*A<=repmat(b,size(R,1),1),2)==size(A,2),:);
        if size(R,1)>0
            X((n+1):min(N,(n+size(R,1))),:) = R(1:min(N-n,size(R,1)),:);
            nar = nar + min(N,(n+size(R,1))) - n;
        end
        n = n + size(R,1); trials = trials + s;
        rho = n/trials;
        if rho>0
            s = min([maxSample,ceil((N-n)/rho),10*s]);
        else
            s = min([maxSample,10*s]);
        end
        passes = passes + 1;
    end
end

%
% Approach 2: Gibbs sampler of Christian Robert
%
if nar < N
    % choose starting point
    if nar>0
        x = X(nar,:);
    else
        x = mu;
    end
    % set up inverse Sigma
    SigmaInv = inv(sigma);
    n = nar;
    while n<N
        % choose p new components
        for i = 1:p
            % Sigmai_i is the (p ? 1) vector derived from the i-th column of ��
            % by removing the i-th row term.
            Sigmai_i = sigma([1:(i-1) (i+1):p],i);
            % Sigma_i_iInv is the inverse of the (p?1)��(p?1) matrix
            % derived from �� = (��ij ) by eliminating its i-th row and
            % its i-th column 
            Sigma_i_iInv = SigmaInv([1:(i-1) (i+1):p],[1:(i-1) (i+1):p]) - ...
                SigmaInv([1:(i-1) (i+1):p],i)*SigmaInv([1:(i-1) (i+1):p],i)' ...
                / SigmaInv(i,i);
            % x_i is the (p-1) vector of components not being updated
            % at this iteration. /// mu_i
            x_i = x([1:(i-1) (i+1):p]);
            mu_i = mu([1:(i-1) (i+1):p]);
            % mui is E(xi|x_i)
            %        mui = mu(i) + Sigmai_i' * Sigma_i_iInv * (x_i - mu_i);
            mui = mu(i) + Sigmai_i' * Sigma_i_iInv * (x_i' - mu_i');
            s2i = sigma(i,i) - Sigmai_i'*Sigma_i_iInv*Sigmai_i;
            % Find points where the line with the (p-1) components x_i
            % fixed intersects the bounding polytope.
            % A_i is the (p-1) x m matrix derived from A by removing
            % the i-th row.
            A_i = A([1:(i-1) (i+1):p],:);
            % Ai is the i-th row of A
            Ai = A(i,:);
            c = (b-x_i*A_i)./Ai;
            lb = max(c(Ai<0));
            if isempty(lb), lb=-Inf; end
            ub = min(c(Ai>0));
            if isempty(ub), ub=Inf; end
            
            %% test feasibility of using TruncatedGaussian
            % lbsig = (lb-mui)/sqrt(s2i);
            % ubsig = (ub-mui)/sqrt(s2i);
            % if lbsig > 10 || ubsig < -10
            %     fprintf('YOWZA! %.2f %.2f %.2f %.2f\n', lbsig, ubsig, ...
            %             mui, sqrt(s2i));
            % end
            %%
            
            % now draw from the 1-d normal truncated to [lb, ub]
            x(i) = mui+TruncatedGaussian(-sqrt(s2i),[lb ub]-mui);
        end
        n = n + 1;
        X(n,:) = x;
        ngibbs = ngibbs+1;
    end
end