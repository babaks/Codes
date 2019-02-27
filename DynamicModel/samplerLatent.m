function output = samplerLatent( y, sigma, Tr, Ind )
[D,n,N] = size(y);
output = zeros(D,n,N);
parfor jj =1:N
    y_temp = zeros(D,n);
	y_temp(:,:) = y(:,:,jj);
	output_jj = zeros(D,n); 
    for ii = 1:n
        Ind_ii = Ind(:,ii);
        I_indep = find( Ind_ii==0 );
        I_corr = find( Ind_ii==1 );
        n_corr = length(I_corr);
        n_indep = length(I_indep);
        if n_corr <= 1
            for kk = 1:D
                if y_temp(kk,ii) == 1
                    output_jj(kk,ii) = rmvnrnd(0,sigma(kk,kk),1,1,Tr(ii,kk),0.5);
                else
                    output_jj(kk,ii) = rmvnrnd(0,sigma(kk,kk),1,-1,-Tr(ii,kk),0.5);
                end
            end
        else     
            mu_corr=zeros(n_corr,1);
            sigma_corr = sigma(I_corr,I_corr);
            if n_indep > 0
                for kk = 1:n_indep
                    if y_temp(I_indep(kk),ii) == 1
                        output_jj(I_indep(kk),ii) = rmvnrnd(0,sigma(I_indep(kk),I_indep(kk)),1,1,Tr(ii,I_indep(kk)),0.5);
                    else
                        output_jj(I_indep(kk),ii) = rmvnrnd(0,sigma(I_indep(kk),I_indep(kk)),1,-1,-Tr(ii,I_indep(kk)),0.5);
                    end
                end                
                A = diag( y_temp(I_corr,ii) .* ones(n_corr,1) - (1-y_temp(I_corr,ii)) .* ones(n_corr,1) );
                b = y_temp(I_corr,ii) .* Tr(ii,I_corr)' - (1-y_temp(I_corr,ii)) .* Tr(ii,I_corr)';
                [X] = rmvnrnd(mu_corr,sigma_corr,1,A,b,0.5);
                output_jj(I_corr,ii) = X;
            else
                A = diag( y_temp(I_corr,ii) .* ones(n_corr,1) - (1-y_temp(I_corr,ii)) .* ones(n_corr,1) );
                b = y_temp(I_corr,ii) .* Tr(ii,I_corr)' - (1-y_temp(I_corr,ii)) .* Tr(ii,I_corr)';
                [X] = rmvnrnd(mu_corr,sigma_corr,1,A,b,0.5);
                output_jj(I_corr,ii) = X;
            end
        end      
    end
    output(:,:,jj) = output_jj;
end

end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    