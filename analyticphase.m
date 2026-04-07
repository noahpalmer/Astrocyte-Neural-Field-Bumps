function phase = analyticphase(beta_vals,D_vals,gamma_vals,theta,choice)
    
    % This function computes a matrix whose entries are estimates of the
    % derivative of the Evans function at \lambda=0. The boundary between
    % regions where stationary bumps are unstable and stable is found by
    % estimating the location of the zero contour.

    % Inputs:
        
        % beta_vals: synaptic depltion values, a vector
        % D_vals: astrocytic diffusion values, a vector or a scalar (for D
        % \to\infty) limit D can be set to any value.
        % gamma_vals: synaptic replenishment values, a vector or a scalar 
        % theta: activity threshold, a scalar
        % choice: set to 1 for (D,\beta) or (\gamma,\beta) space with finite D, set to 2 for (\gamma,\beta) space for large
        % diffusion limit.

    % Note: one of gamma_vals or D_vals must be a scalar.

    % Outputs:

        % phase: matrix of approximations to E' at \lambda=0.

    if length(gamma_vals) > 1 && length(D_vals) > 1
        error('The function "analyticphase" can only be used to create data for one phase diagram at a time. As a consequence, one of gamma_vals or D_vals must be a scalar.');
    end

    gammaN = length(gamma_vals);
    DN = length(D_vals);
    betaN = length(beta_vals);
    

    % step size for centered difference approximation to E'.
    h = 1e-12;

    if length(D_vals) > 1
        gamma = gamma_vals;
        phase = zeros(DN,betaN);
            for i = 1:DN
                D = D_vals(i);
                for j = 1:betaN
                    beta = beta_vals(j);
            
                    f = @(delta) czero(delta,beta,gamma).*sin(2*delta)-theta;
                    a = fzero(f,1.5);
                    
                    % Auxiliary parameters
                    kappa = a/pi;
                    c0 = czero(a,beta,gamma);
                    A0 = kappa*(1-c0);
                    mu = (beta+gamma*A0)/(2*gamma*A0*(sin(a))^2);
                    Qp = 1;
                    Qm = c0;
            
                    % Centered difference approximation to E'
                    Dplus  = detcomp(h,D,beta,gamma,A0,c0,mu,a,Qp,Qm,1);
                    Dminus = detcomp(-h,D,beta,gamma,A0,c0,mu,a,Qp,Qm,1);
                    
                    Dprime = (Dplus-Dminus)/(2*h);
                    
                    phase(j,i) = real(Dprime);
                end 
            end

    % Handles both the finite and infinite diffusion cases depending on
    % choice value.
    else
            D = D_vals;
            phase = zeros(gammaN,betaN);
            for i = 1:gammaN
                gamma = gamma_vals(i);
        
                for j = 1:betaN
                    beta = beta_vals(j);
    
                    f = @(delta) czero(delta,beta,gamma).*sin(2*delta)-theta;
                    a = fzero(f,1.5);
                    
                    % Auxiliary parameters
                    kappa = a/pi;
                    c0 = czero(a,beta,gamma);
                    A0 = kappa*(1-c0);
                    mu = (beta+gamma*A0)/(2*gamma*A0*(sin(a))^2);
                    Qp = 1;
                    Qm = c0;
    
                    % Centered difference approximation to E'
                    gammaplus  = detcomp(h,D,beta,gamma,A0,c0,mu,a,Qp,Qm,choice);
                    gammaminus = detcomp(-h,D,beta,gamma,A0,c0,mu,a,Qp,Qm,choice);
            
                    gammaprime = (gammaplus-gammaminus)/(2*h);
            
                    phase(j,i) = real(gammaprime);
                end  
            end
    end
end

function c0 = czero(delta,beta,gamma)  
    c0 = (beta+2*(gamma*delta/pi)-sqrt(beta^2+4*beta*gamma*delta/pi))./(2*(gamma*delta/pi));
end
