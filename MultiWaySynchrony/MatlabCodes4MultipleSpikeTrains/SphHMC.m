%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to sample latent variables using Spherical HMC
% Inputs: current_q is the current state of copula model parameters, UUU is the function to calculate log-pmf and derivative of log-pmf, eps and L are the step size and step number
% Outputs: New state of copula model parameters and indicator of acceptance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [q, acpt] = SphHMC(current_q,UUU,eps,L)

q = current_q;
v = randn(length(q),1); % sample velocity
v = v - q*(q'*v); % project to tangent plane
current_E = UUU(q,0) + (v'*v)/2;% current energy


randL=ceil(rand*L); % Number of leap frog

%% forward half step of velocity
v = v - transpose(eps/2*UUU(q,1));

for i=1:randL
% forward half step of velocity
%     v = v - transpose(eps/2*UUU(q,1));
%     if(abs(q'*v)>1e-6)
%         v = v - q*(q'*v); % calibrate direction possibly deviated by error accumulation
%    end
    
    %% full step evolution on sphere along great circle
    q0 = q; v_nom = sqrt(v'*v);
    q = q0*cos(v_nom*eps) + v/v_nom*sin(v_nom*eps);
    v = -q0*v_nom*sin(v_nom*eps) + v*cos(v_nom*eps);
    
    %% Full step of velocity
    if i~=randL
        v = v - transpose(eps*UUU(q,1));
    end    
end

%% Forward half step of velocity
 v = v - transpose(eps/2*UUU(q,1));

proposed_E = UUU(q,0) + (v'*v)/2; % new energy

%% log of Metropolis ratio
logRatio = -proposed_E + current_E;
if (isfinite(logRatio) & (log(rand) < min([0,logRatio])))
    acpt = 1;
    q=q;
else
    q = current_q;
    acpt = 0;
end

end