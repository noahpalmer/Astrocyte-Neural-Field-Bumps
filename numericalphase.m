function [velocity, drift] = numericalphase(beta_vals,D_vals,gamma_vals,theta,N,T)
    % This function computes a matrix whose entries are estimates of the
    % derivative of the Evans function at \lambda=0. The boundary between
    % regions where stationary bumps are unstable and stable is found by
    % estimating the location of the zero contour.

    % Inputs:
        
        % beta_vals: synaptic depletion values, a vector
        % gamma_vals: synaptic replenishment, a vector or a scalar
        % D_vals: astrocytic diffusion values, a vector or a scalar 
        % theta: activity threshold, a scalar
        % N: number of spatial grid points.
        % T: final time t for the simulations.

    % Note: one of gamma_vals or D_vals must be a scalar.

    % Outputs:

        % velocity: matrix containing approximations to the velocity of the
        % center of mass of the bump at time T.
        % drift: matrix containing approximations to the cumulative drift
        % of the center of mass of the bump at time T.

    if length(gamma_vals) > 1 && length(D_vals) > 1
        error('The function "numericalphase" can only be used to create data for one phase diagram at a time. As a consequence, one of gamma_vals or D_vals must be a scalar.');
    end
    
    starttime = tic;
   
    dx = 2*pi/N;
    x = linspace(-pi,pi-dx,N)';

    dt = 0.01;
    nt = round(T/dt)+1;

    e = ones(N,1);
    D2 = spdiags([e -2*e e],-1:1,N,N);
    D2(1,end) = 1;
    D2(end,1) = 1;
    D2 = D2/dx^2;
    IA = speye(N);
    
    cosx = cos(x);
    sinx = sin(x);
    cx = dx*cosx'; 
    sx = dx*sinx';

    if length(D_vals) > 1
        velocity = zeros(length(beta_vals),length(D_vals));
        drift = zeros(length(beta_vals),length(D_vals));
        gamma = gamma_vals;
        totaliterations = length(beta_vals)*length(D_vals);
        iter = 0;
        for i = 1:length(D_vals)
            D = D_vals(i);

            LA = IA-(dt*D)*D2;
            LAdec = decomposition(LA,'lu');
        
            Delta_prev = 1.5; % Initial guess is for larger bump
        
            for j = 1:length(beta_vals)
        
                beta = beta_vals(j);
        
                % Solve for $\Delta$
                f = @(delta) czero(delta,beta,gamma).*sin(2*delta)-theta;
                Delta = fzero(f,Delta_prev);
                Delta_prev = Delta;
        
                kappa = Delta/pi;
                c0 = czero(Delta,beta,gamma);
                A0 = kappa*(1-c0);
                epsilon = 0.05*(2*c0*sin(Delta));  %Perturbation size epsilon taken to be 0.05 of maximum height stationary solution
        
                % Initial conditions
                U = (2*c0*sin(Delta))*cosx+epsilon*sin(x);
                A = A0;
                Q = c0*((x>-Delta)&(x<Delta))+((x>Delta)|(x<-Delta));
        
                % Velocity is determined by tracking the center of mass between
                % time steps using a circular mean.
                % So we initialize the center of mass for the IC:
                mass = trapz(x,U);
                C0 = (cx*U)/mass;
                S0 = (sx*U)/mass;
                com_start = atan2(S0,C0);
                com_old = com_start;
                com_total = 0;

        
                for k = 1:nt-1
        
                    H = (U > theta);
                    QHu = Q.*H;
        
                    fcos = cx*QHu;
                    fsin = sx*QHu;
                    Unew = (1-dt)*U+dt*(fcos*cosx+fsin*sinx);
        
                    % Update center of mass and compute speed.
                    
                    C = (cx*Unew)/mass;
                    S = (sx*Unew)/mass;
                    com_new = atan2(S,C);

                    dcom = com_new-com_old;
                    if dcom > pi
                        dcom = dcom-2*pi;
                    elseif dcom < -pi
                        dcom = dcom+2*pi;
                    end

                    com_total = com_total+dcom;
                    com_old = com_new;

                    c = gamma*A;
                    lam = c+beta*H;
        
                    E = exp(-(dt)*lam);
                    Qnew = E.*Q+(1-E).*(c./lam);
        
                    Rexp = (dt)*(beta*QHu-gamma*A.*(1-Q));
                    Anew = LAdec\(A+Rexp);
        
                    U = Unew;
                    A = Anew;
                    Q = Qnew;
                end
                velocity(j,i) = abs(dcom)/dt;
                drift(j,i) = com_total;

                iter = iter + 1;
                elapsedtime = toc(starttime);
                percent = (iter/totaliterations)*100;
                totaltime = elapsedtime/iter*totaliterations;
                remainingtime = totaltime-elapsedtime;

                fprintf('\rProgress: %6.2f%% | Elapsed: %6.1fs | Estimated Remaining Time: %6.1fs',percent, elapsedtime, remainingtime);
            end
        end
    else
        velocity = zeros(length(beta_vals),length(gamma_vals));
        drift = zeros(length(beta_vals),length(gamma_vals));
        D = D_vals;
        totaliterations = length(beta_vals)*length(gamma_vals);

        LA = IA-(dt*D)*D2;
        LAdec = decomposition(LA,'lu');

        iter = 0;

        for i = 1:length(gamma_vals)

                gamma = gamma_vals(i);
            
                Delta_prev = 1.5; % Initial guess is for larger bump 
            
                for j = 1:length(beta_vals)
                   
                    beta = beta_vals(j);
            
                    % Solve for $\Delta$
                    f = @(delta) czero(delta,beta,gamma).*sin(2*delta)-theta;
                    Delta = fzero(f,Delta_prev);
                    Delta_prev = Delta;
            
                    kappa = Delta/pi;
                    c0 = czero(Delta,beta,gamma);
                    A0 = kappa*(1-c0);
                    epsilon = 0.05*(2*c0*sin(Delta));  %Perturbation size epsilon taken to be 0.05 of maximum height stationary solution
            
                    % Initial conditions
                    U = (2*c0*sin(Delta))*cosx+epsilon*sin(x);
                    A = A0;
                    Q = c0*((x>-Delta)&(x<Delta))+((x>Delta)|(x<-Delta));
            
                    % Velocity is determined by tracking the center of mass between
                    % time steps using a circular mean.
                    % So we initialize the center of mass for the IC:
                    mass = trapz(x,U);
                    C0 = (cx*U)/mass;
                    S0 = (sx*U)/mass;
                    com_start = atan2(S0,C0);
                    com_old = com_start;
                    com_total = 0;
            
                    for k = 1:nt-1
            
                        H = (U > theta);
                        QHu = Q.*H;
            
                        fcos = cx*QHu;
                        fsin = sx*QHu;
                        Unew = (1-dt)*U+dt*(fcos*cosx+fsin*sinx);
            
                        % Update center of mass and compute speed.
                        C = (cx*Unew)/mass;
                        S = (sx*Unew)/mass;
                        com_new = atan2(S,C);
            
                        dcom = com_new-com_old;
                        if dcom > pi
                            dcom = dcom-2*pi;
                        elseif dcom < -pi
                            dcom = dcom+2*pi;
                        end
            
                        com_total = com_total+dcom;
                        com_old = com_new;
            
                        c = gamma*A;
                        lam = c+beta*H;
            
                        E = exp(-(dt)*lam);
                        Qnew = E.*Q+(1-E).*(c./lam);
            
                        Rexp = (dt)*(beta*QHu-gamma*A.*(1-Q));
                        Anew = LAdec\(A+Rexp);
            
                        U = Unew;
                        A = Anew;
                        Q = Qnew;
                    end
                    velocity(j,i) = abs(dcom)/dt;
                    drift(j,i) = com_total;

                    iter = iter + 1;
                    elapsedtime = toc(starttime);
                    percent = (iter/totaliterations)*100;
                    totaltime = elapsedtime/iter*totaliterations;
                    remainingtime = totaltime-elapsedtime;

                    fprintf('\rProgress: %6.2f%% | Elapsed: %6.1fs | Estimated Remaining Time: %6.1fs',percent, elapsedtime, remainingtime);
                end
         end
    end
    fprintf('\nTime elapsed: %.0f seconds\n',toc(starttime));
end


function c0 = czero(delta,beta,gamma)  
    c0 = (beta+2*(gamma*delta/pi)-sqrt(beta^2+4*beta*gamma*delta/pi))./(2*(gamma*delta/pi));
end

