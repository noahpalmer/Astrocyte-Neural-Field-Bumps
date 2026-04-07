function  d = detcomp(lambda,D,beta,gamma,A0,c0,mu,a,Qp,Qm,choice)
% This function computes the Evans function for the given parameters and
% lambda choice. It is used in "analyticphase.m".

    % Inputs:
        % beta: synaptic depletion
        % D: astrocytic diffusion
        % gamma: synaptic replenishment
        % A0: amplitude of astrocyte resource pool for stationary bump
        % c0: amplitude of synaptic resource pool in the active region for
        % stationary bump
        % Qp: value of Q_+, case dependent
        % Qm: value of Q_-, case dependent
        % mu: 1/|U'(Delta)|
        % a: bump width
        % choice: set to 1 finite diffusion case, set to 2 for large
        % diffusion limit

    % Outputs:

        % d: value of the Evans function for the given parameters

    E = [ (Qp+Qm)*cos(a)^2, (Qp-Qm)*cos(a)*sin(a); 
          (Qp-Qm)*cos(a)*sin(a), (Qp+Qm)*sin(a)^2 ];
    if choice == 2 
        B = Bcompinf(lambda,beta,gamma,c0,A0,mu,a,Qp,Qm); % Large diffusion limit
        M = (lambda+1)*eye(2)-mu*E-B;
        d = det(M);
    else
        B = Bcompgen(lambda,D,beta,gamma,A0,c0,mu,a,Qp,Qm); % Finite diffusion
        M = (lambda+1)*eye(2)-mu*E-B;
        d = det(M);
    end
end


function [B,gplus,gminus] = Bcompgen(lambda,D,beta,gamma,A0,c0,mu,a,Qp,Qm)

    c1 = lambda/D;
    c2 = lambda*(lambda+gamma*A0+beta+gamma*(1-c0))/(D*(lambda+gamma*A0+beta));
    c3 = lambda/D;

    s1 = sqrt(c1 + 0i);
    s2 = sqrt(c2 + 0i);
    s3 = sqrt(c3 + 0i);

    Ap = zeros(6,6);
    bp = zeros(6,1);
    Am = zeros(6,6);
    bm = zeros(6,1);

    Ap(1,:) = [exp(-s1*pi), exp(s1*pi), 0,0, -exp(s3*pi), -exp(-s3*pi)];
    Ap(2,:) = [s1*exp(-s1*pi), -s1*exp(s1*pi), 0,0, -s3*exp(s3*pi), s3*exp(-s3*pi)];
    Ap(3,:) = [exp(-s1*a), exp(s1*a), -exp(-s2*a), -exp(s2*a), 0,0];
    Ap(4,:) = [s1*exp(-s1*a), -s1*exp(s1*a),-s2*exp(-s2*a),s2*exp(s2*a),0,0];
    Ap(5,:) = [0, 0, exp(s2*a), exp(-s2*a), -exp(s3*a), -exp(-s3*a)];
    Ap(6,:) = [0,0, -s2*exp(s2*a), s2*exp(-s2*a), s3*exp(s3*a), -s3*exp(-s3*a)];
    bp(6) = -((mu*beta*lambda))/(D*(lambda+gamma*A0+beta));

    coeffp = Ap\bp;

    Am(1,:) = [exp(-s1*pi), exp(s1*pi), 0,0, -exp(s3*pi), -exp(-s3*pi)];
    Am(2,:) = [s1*exp(-s1*pi), -s1*exp(s1*pi), 0,0, -s3*exp(s3*pi), s3*exp(-s3*pi)];
    Am(3,:) = [exp(-s1*a), exp(s1*a), -exp(-s2*a), -exp(s2*a), 0,0];
    Am(4,:) = [-s1*exp(-s1*a), s1*exp(s1*a), s2*exp(-s2*a), -s2*exp(s2*a),0,0];
    Am(5,:) = [0, 0, exp(s2*a), exp(-s2*a), -exp(s3*a), -exp(-s3*a)];
    Am(6,:) = [0,0, s2*exp(s2*a), -s2*exp(-s2*a), -s3*exp(s3*a), s3*exp(-s3*a)];
    bm(4) = -((mu*beta*lambda))/(D*(lambda+gamma*A0));

    coeffm = Am\bm;
    gplus  = @(y) coeffp(3)*exp(s2*y)+coeffp(4)*exp(-s2*y);
    gminus = @(y) coeffm(3)*exp(s2*y)+coeffm(4)*exp(-s2*y);


    fB11 = @(y) cos(y).*(cos(a).*(Qp.*gplus(y)+Qm.*gminus(y)));
    fB12 = @(y) cos(y).*(sin(a).*(Qp.*gplus(y)-Qm.*gminus(y)));
    fB21 = @(y) sin(y).*(cos(a).*(Qp.*gplus(y)+Qm.*gminus(y)));
    fB22 = @(y) sin(y).*(sin(a).*(Qp.*gplus(y)-Qm.*gminus(y)));

    B11 = integral(fB11,-a,a,'ArrayValued',true);
    B12 = integral(fB12,-a,a,'ArrayValued',true);
    B21 = integral(fB21,-a,a,'ArrayValued',true);
    B22 = integral(fB22,-a,a,'ArrayValued',true);

    C = (gamma*(1-c0))/(lambda+beta+gamma*A0)*[B11 B12; B21 B22]; 

    A = -((mu*beta*Qp)/(lambda+beta+gamma*A0))*[(cos(a))^2 cos(a)*sin(a); cos(a)*sin(a) (sin(a))^2];
    B = C+A;
end


function B = Bcompinf(lambda,beta,gamma,c0,A0,mu,a,Qp,Qm)

    K1 = (gamma*(1-c0)*Qp)/(pi*(lambda+gamma*A0+beta)+a*gamma*(1-c0));
    K2 = (gamma*(1-c0)*Qm*(lambda+gamma*A0+beta))/((lambda+gamma*A0)*(pi*(lambda+gamma*A0+beta)+a*gamma*(1-c0)));
    
    B11 = -((mu*beta)/(lambda+gamma*A0+beta))*((cos(a))^2+(K1+K2)*sin(a)*cos(a));
    B12 = -((mu*beta)/(lambda+gamma*A0+beta))*(cos(a)*sin(a)+(K1-K2)*(sin(a))^2);
    B21 = -((mu*beta)/(lambda+gamma*A0+beta))*(cos(a)*sin(a));
    B22 = -((mu*beta)/(lambda+gamma*A0+beta))*(sin(a))^2;

    B = [B11 B12; B21 B22];
end


