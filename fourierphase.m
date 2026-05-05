function stablevals = fourierphase(beta_vals,D_vals,gamma_vals,theta,tol)

    % This function computes a matrix whose entries are estimates of the
    % derivative of the Evans function at \lambda=0. The boundary between
    % regions where stationary bumps are unstable and stable is found by
    % estimating the location of the zero contour.

    % Inputs:
        
        % beta_vals: synaptic depltion values, a vector
        % D_vals: astrocytic diffusion values, a vector or a scalar 
        % gamma_vals: synaptic replenishment values, a vector or a scalar 
        % theta: activity threshold, a scalar
        % tol: to classify whether a numerically computed eigenvalue is
        % positive or not.
        % initial_guess: 
        % choice: set to 1 for (D,\beta) space, set to 2 for (\gamma,\beta)
        % space with finite D

    % Note: one of gamma_vals or D_vals must be a scalar.

    % Outputs:
        
        % stablevals: a matrix where 1s denote at least one eigenvalue is positive
        % for those parameters and 0s denote all eigenvalues are less
        % than or equal to 0.

    if length(gamma_vals) > 1 && length(D_vals) > 1
        error('The function "fourierphase" can only be used to create data for one phase diagram at a time. As a consequence, one of gamma_vals or D_vals must be a scalar.');
    end


    M = [1, 0, 0, 0, 0, 0;
        0, 1, 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, 1, 0, 0;
        0, 0, -1, 0, 1, 0;
        0, 0, 0, -1, 0, 1];
    initial_guess = 1.5; % initial guess of bump width, set for larger bump
    
    
    if length(D_vals) > 1
        gamma = gamma_vals;
        stablevals = false(length(beta_vals),length(D_vals));
        Deltavals = zeros(1,length(beta_vals));

        for i = 1:length(beta_vals)
            beta = beta_vals(i);
            f = @(delta) czero(delta,beta,gamma).*sin(2*delta)-theta;
            
            Deltavals(i) = fzero(f,initial_guess);
            Delta = Deltavals(i);
            kappa = Delta/pi;
            initial_guess = Delta;
        
            % Compute amplitude of Q in the active region, the value of A and
            % 1/U'(\Delta).
            c0 = (beta+2*gamma*kappa-sqrt(beta^2+4*beta*gamma*kappa))/(2*gamma*kappa);
            A0 = kappa*(1-c0);
            mu = (beta+gamma*A0)/(2*gamma*A0*(sin(Delta))^2);

            % Perturbation specific coefficients
            f1 = mu*(1+c0)*(cos(Delta))^2;
            f2 = mu*(1-c0)*cos(Delta)*sin(Delta);
            f3 = mu*(1+c0)*(sin(Delta))^2;
            f4 = mu*beta*(cos(Delta))^2;
            f5 = mu*beta*cos(Delta)*sin(Delta);
            f6 = mu*beta*(sin(Delta))^2;
    
            for k = 1:length(D_vals)
                D = D_vals(k);
            
                A = [f1-1, f2, 1,  0,  0, 0;
                    f2, f3-1,  0,  1,  0, 0;
                    -f4, -f5, -(gamma*A0+beta), 0, gamma*(1-c0), 0;
                    -f5, -f6, 0, -(gamma*A0+beta), 0, gamma*(1-c0);
                    0, 0, 0, 0, -D, 0;
                    0, 0, 0, 0, 0, -D];
    
                Z = M*A;
                eigvals = eig(Z);
                [V,~] = eig(Z);
                
                % Test that the eigenvectors satisfy the right shift condition.
                % Keep eigenvalues that meet this assumption.
    
                valid_eigvals = [];
                for j = 1:length(eigvals)
                    v = V(:,j); 
                    test1 = (v(1)*cos(Delta)+v(2)*sin(Delta) > 0) && (v(1)*cos(Delta)-v(2)*sin(Delta) < 0);
                    test2 = (-v(1)*cos(Delta)-v(2)*sin(Delta) > 0) && (-v(1)*cos(Delta)+v(2)*sin(Delta) < 0);
                    if (test1 || test2)
                        valid_eigvals(end+1) = eigvals(j);  
                    end
                end
                if any(real(valid_eigvals) > tol)
                    stablevals(i,k) = true;
                end
            end
        end
    else
        D = D_vals;
        stablevals = false(length(beta_vals),length(gamma_vals));
        for i = 1:length(beta_vals)
            beta = beta_vals(i);
        
            for k = 1:length(gamma_vals)
                gamma = gamma_vals(k);
        
        
                f = @(delta) czero(delta,beta,gamma).*sin(2*delta)-theta;
                Delta = fzero(f,initial_guess);
                initial_guess = Delta;
        
                kappa = Delta/pi;
        
                c0 = (beta+2*gamma*kappa-sqrt(beta^2+4*beta*gamma*kappa))/(2*gamma*kappa);
                A0 = kappa*(1-c0);
                mu = (beta+gamma*A0)/(2*gamma*A0*(sin(Delta))^2);
        
                % Perturbation specific coefficients
                f1 = mu*(1+c0)*(cos(Delta))^2;
                f2 = mu*(1-c0)*cos(Delta)*sin(Delta);
                f3 = mu*(1+c0)*(sin(Delta))^2;
                f4 = mu*beta*(cos(Delta))^2;
                f5 = mu*beta*cos(Delta)*sin(Delta);
                f6 = mu*beta*(sin(Delta))^2;
        
                A = [f1-1, f2, 1,  0,  0, 0;
                        f2, f3-1,  0,  1,  0, 0;
                        -f4, -f5, -(gamma*A0+beta), 0, gamma*(1-c0), 0;
                        -f5, -f6, 0, -(gamma*A0+beta), 0, gamma*(1-c0);
                        0, 0, 0, 0, -D, 0;
                        0, 0, 0, 0, 0, -D];
                Z = M*A;
        
                [V,E] = eig(Z);
                eigvals = diag(E);   
        
                % Test that the eigenvectors satisfy the right shift condition.
                % Keep eigenvalues that meet this assumption.
                valid_eigvals = [];
                for j = 1:length(eigvals)
                    v = V(:,j);
        
                    test1 = (v(1)*cos(Delta)+v(2)*sin(Delta) > 0) && (v(1)*cos(Delta)-v(2)*sin(Delta) < 0);
                    test2 = (-v(1)*cos(Delta)-v(2)*sin(Delta) > 0) && (-v(1)*cos(Delta) + v(2)*sin(Delta) < 0);
        
                    if test1 || test2
                        valid_eigvals(end+1) = eigvals(j); 
                    end
                end
        
                if any(real(valid_eigvals) > tol)
                    stablevals(i,k) = true;
                end
            end
        end
    end
end

function c0 = czero(delta,beta,gamma)  
    c0 = (beta+2*(gamma*delta/pi)-sqrt(beta^2+4*beta*gamma*delta/pi))./(2*(gamma*delta/pi));
end
