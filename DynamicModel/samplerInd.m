function [output1,output2] = samplerInd ( uu, Sigma, ind_neuron1, ind_neuron2, Ind, pt )
[D,n] = size(Ind);
output1 = zeros(n,1);
output2 = zeros(n,1);

parfor ii = 1:n
    Ind11 = Ind(:,ii);
    Ind11(ind_neuron1) = 1;
    Ind11(ind_neuron2) = 1;
    loglike11 = loglikeSigma( uu(:,ii,:), Sigma, Ind11 );
    prior11 = log(pt(ind_neuron1,ii)+0.00000001) + log(pt(ind_neuron2,ii)+0.00000001);

    Ind10 = Ind(:,ii);
    Ind10(ind_neuron1) = 1;
    Ind10(ind_neuron2) = 0;
    loglike10 = loglikeSigma( uu(:,ii,:), Sigma, Ind10 );
    prior10 = log(pt(ind_neuron1,ii)+0.00000001) + log(1-pt(ind_neuron2,ii)+0.00000001);

    Ind01 = Ind(:,ii);
    Ind01(ind_neuron1) = 0;
    Ind01(ind_neuron2) = 1;
    loglike01 = loglikeSigma( uu(:,ii,:), Sigma, Ind01 );
    prior01 = log(1-pt(ind_neuron1,ii)+0.00000001) + log(pt(ind_neuron2,ii)+0.00000001);

    Ind00 = Ind(:,ii);
    Ind00(ind_neuron1) = 0;
    Ind00(ind_neuron2) = 0;
    loglike00 = loglikeSigma( uu(:,ii,:), Sigma, Ind00 );
    prior00 = log(1-pt(ind_neuron1,ii)+0.00000001) + log(1-pt(ind_neuron2,ii)+0.00000001);
    

    p11 = 1/(1+exp(loglike10+prior10-loglike11-prior11)+exp(loglike01+prior01-loglike11-prior11)+exp(loglike00+prior00-loglike11-prior11));
    p10 = 1/(1+exp(loglike11+prior11-loglike10-prior10)+exp(loglike01+prior01-loglike10-prior10)+exp(loglike00+prior00-loglike10-prior10));
	p01 = 1/(1+exp(loglike11+prior11-loglike01-prior01)+exp(loglike10+prior10-loglike01-prior01)+exp(loglike00+prior00-loglike01-prior01));
    p00 = 1/(1+exp(loglike11+prior11-loglike00-prior00)+exp(loglike01+prior01-loglike00-prior00)+exp(loglike10+prior10-loglike00-prior00));
    output1(ii) = binornd(1,p11+p10);
    if output1(ii) == 1
        output2(ii) = binornd(1,p11/(p11+p10));
        if isnan(output2(ii))
            disp('oops!');
        end
    else
        output2(ii) = binornd(1,p01/(p00+p01));
        if isnan(output2(ii))
            disp('oops!');
        end
    end
    
    
% %     disp( [p00,p01,p10,p11] );
%     mn_samp = mnrnd(1, [p00,p01,p10,p11] );    
%     output1(ii,:) = mn_samp(3)==1 | mn_samp(4)==1;
%     output2(ii,:) = mn_samp(2)==1 | mn_samp(4)==1;
end

end
