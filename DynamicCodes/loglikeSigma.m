function output = loglikeSigma( uu, sigma, Ind )

[D,n,N] = size(uu);
loglike = zeros( n, 1 );


for ii = 1:n
    sigma_ii = zeros(D,D);
    for jj = 1:(D-1)
        for kk = (jj+1):D
            sigma_jk = sigma(jj,kk)*Ind(jj,ii)*Ind(kk,ii);
            sigma_ii(jj,kk) = sigma_jk;
            sigma_ii(kk,jj) = sigma_jk;
        end
    end
    for jj =1:D
        sigma_ii(jj,jj) = sigma(jj,jj);
    end
    cholSigma = chol( sigma_ii ); % Cholesky decomposition
    logdetSigma = 2 * sum( log(diag(cholSigma)) ); % Logarithm of the determinant of the covariance matrix
    invChol = inv( cholSigma ); % Inverse of LL


    loglikeii = 0;
    uu_temp = zeros(D,N);
    uu_temp(:,:) = uu(:,ii,:);
    for jj = 1:N
      	invLLuu = transpose(invChol) * uu_temp(:,jj);
        loglikeii = loglikeii - 0.5*logdetSigma - 0.5*transpose(invLLuu)*invLLuu;
    end
    loglike(ii) = loglikeii;
end

output = sum(loglike);

end
