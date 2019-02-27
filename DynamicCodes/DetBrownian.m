function [detSigma] = DetBrownian( t, n )
detSigma =exp( sum(log(t(1:n)-[0,t(1:(n-1))])) );
end
