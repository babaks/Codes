%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate pmf and derivative of pmf for Spike Train Data
% Inputs: yy is the data  and pp is the firing probabilities
% Outputs: pmf1smdDATA, which is not pmf but will be used to calculate pmf, and derivative of pmf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pmf1smdDATA,dpmfDATA] = pmfdpmfDATA(yy,pp)

[N,n,D] = size(yy); % Size of data
pmf1smdtemp = zeros(N,n); % Initialization
dpmftemp = zeros(N,n,D*(D-1)/2); % Initialization for dpmf

parfor ii=1:N % Use parfor to calculate pmf1smdtemp and dpmftemp in parallel
    for jj=1:n
        ytemp = zeros(1,D);
        ptemp = zeros(1,D);
        ytemp(:) = yy(ii,jj,:);
        ptemp(:) = pp(:,jj);
        [output] = pmfdpmf(ytemp,ptemp); 
        pmf1smdtemp(ii,jj) = output(1);
        dpmftemp(ii,jj,:) = output(2:end);
    end
end

pmf1smdDATA = zeros(N*n,1);
dpmfDATA = zeros(N*n,D*(D-1)/2);
for ii = 1:N
    for jj = 1:n
        pmf1smdDATA((ii-1)*n+jj) = pmf1smdtemp(ii,jj);
        dpmfDATA((ii-1)*n+jj,:) = dpmftemp(ii,jj,:);
    end
end

end

