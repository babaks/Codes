%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculatate potential and its gradient w.r.t. theta
% Inputs: theta_til is the current state of theta, der is an indicator(0-potential,1-derivative), pmf1smdDATA is used to calculate pmf, dpmfDATA is the derivative of pmf
% Outputs: Potential or its gradient w.r.t. theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[output]=Ut(theta_til,der,pmf1smdDATA,dpmfDATA)

THETA = theta_til(1:end-1);  
BETA = abs(THETA).*THETA; % Map theta to beta
if der==0
    output = U(BETA,0,pmf1smdDATA,dpmfDATA); % Potential
elseif der==1
    dU = U(BETA,1,pmf1smdDATA,dpmfDATA).*2.*abs(THETA'); % D(D-1)/2 vector
    dUt = [dU 0] - theta_til'.*(dU*THETA); % D(D+1)/2 vector, same size as theta, in tangent space
    output = dUt; % Gradient w.r.t. theta
else
    error('wrong choice of der!');
end

end
