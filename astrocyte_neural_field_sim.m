% Finite difference simulation of astrocyte neural field model.
%
% Runs a simulation of the astrocyte-neural field model for specified initial conditions 
% Currently setup to run a simulation with IC near the stationary bump solution
% derived for the model with a right shift perturbation.


% Spatial and temporal gride and spacing.
N = 3000;    
dx = 2*pi/N;    
x = linspace(-pi,pi-dx,N)'; 

dt = 0.01;  
T = 200;    
nt = round(T/dt)+1;


% Model parameters
gamma = 2; % Synaptic replenishment rate
beta = 0.05; % Synaptic depletion rate
tau = 1; % Synaptic depression timescale
D = 0.7; % Astrocytic resource diffusion constant
theta = 0.1; %Neural activity threshold


% Use threshold condition to find the bump half-width

f = @(delta) czero(delta,beta,gamma).*sin(2*delta)-theta;
initial_guess = 1.5; 
Delta = fzero(f, initial_guess);

% Compute astrocyte resource amplitude and synaptic resource amplitude in the active region for
% stationary bump solutions.

kappa = Delta/pi;
c0 = czero(Delta,beta,gamma);
A0 = kappa*(1-c0);
epsilon = 0.3*(2*c0*sin(Delta));


% Initial conditions: stationary bump profile of section 3

U = zeros(N,nt); Q = ones(N,nt); A = U;
U(:,1) = (2*c0*sin(Delta)).*cos(x)+epsilon.*sin(x);
A(:,1) = A0;
Q(:,1) = (c0).*((x>-Delta) & (x<Delta))+1.*((x>Delta)|(x<-Delta));

% Diffusion operator with periodic BCs

e = ones(N,1);
D2 = spdiags([e -2*e e], -1:1, N, N);
D2(1,end) = 1;  
D2(end,1) = 1;  
D2 = D2/dx^2;
IA = speye(N);
LA = IA - (dt*D/tau)*D2;

% Precompute sine and cosine terms times dx needed for integral term in the equation for u(x,t).

cx = dx*cos(x'); sx = dx*sin(x');


for k=1:nt -1

    %u(x,t) update
    QHu = Q(:,k).*(U(:,k) > theta);
    fcos = cx*QHu; fsin = sx*QHu;
    U(:,k+1) = (1-dt)*U(:,k)+dt*(fcos*cos(x)+fsin*sin(x));
    
    % q(x,t) update, first compute integrating factor for q(x,t) 
    c = gamma*A(:,k);
    lam = c+beta*(U(:,k) > theta);
    E = exp(-(dt/tau)*lam);

    Q(:,k+1) = E.*Q(:,k) + (1 - E).*(c./max(lam,1e-12));

    % a(x,t) update, backward euler
    Rexp = (dt/tau)*(beta*QHu-gamma*A(:,k).*(1-Q(:,k)));
    A(:,k+1) = LA\(A(:,k)+Rexp);
end


function c0 = czero(delta, beta, gamma)  
    c0 = (beta+2*(gamma*delta/pi)-sqrt(beta^2+4*beta*gamma*delta/pi))./(2*(gamma*delta/pi));
end


%%


% Code in this section creates Figures 2 and 3.

% Figure 2

% These values are estimates of where bumps transition from stable to unstable for a
% particular theta value (theta=0.03 for beta_cut1 and theta=0.3 for
% beta_cut2) as predicted by the stability analysis of section 4.

beta_cut1 = 0.0255017;   
beta_cut2 = 0.01987; 


% First curve theta=0.03
theta1 = 0.03;    
betas1 = linspace(0.001, beta_cut1, 600);  

Delta_sols1u = zeros(size(betas1));
Dg1u = 1.5;
for k = 1:length(betas1)
    b = betas1(k);
    
    f = @(D) theta1 - ((b + 2*gamma*D/pi - sqrt(b^2 + 4*b*gamma*D/pi)) / (2*gamma*D/pi)) .* sin(2*D);
    Delta0 = Dg1u;
    Delta_sols1u(k) = fsolve(f, Delta0);
end


beta1  = linspace(0.001, 20, 400);
Deltap1 = linspace(0.01, pi/2, 400);  

[B1, D1] = meshgrid(beta1, Deltap1);
c0 = (B1 + 2*gamma*D1/pi - sqrt(B1.^2 + 4*B1*gamma.*D1/pi))./ (2*gamma*D1/pi);
F1 = theta1-c0.*sin(2*D1);


% Second curve theta=0.3
theta2 = 0.3;    
betas2 = linspace(0.001, beta_cut2, 600);  

Delta_sols2u = zeros(size(betas2));
Dg2u = 1.5;
for k = 1:length(betas2)
    b = betas2(k);
    
    f = @(D) theta2 - ((b + 2*gamma*D/pi - sqrt(b^2 + 4*b*gamma*D/pi)) / (2*gamma*D/pi)) .* sin(2*D);
    Delta0 = Dg2u;
    Delta_sols2u(k) = fsolve(f, Delta0);
end




beta2  = linspace(0.001, 5, 400);
Deltap2 = linspace(0.01, pi/2, 400);   

[B2, D2] = meshgrid(beta2, Deltap2);
c0 = (B2 + 2*gamma*D2/pi - sqrt(B2.^2 + 4*B2*gamma.*D2/pi))./ (2*gamma*D2/pi);
F2 = theta2-c0.*sin(2*D2);


% Figure requires running finite difference simulation to first to obtain
% U(:,1), Q(:,1) and A(:,1).

figure('Color','w','Position',[100 100 900 700])
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')

nexttile
plot(x,Q(:,1),'k-.','LineWidth',5); hold on;
plot(x,A(:,1),'k--','LineWidth',5)
plot(x,U(:,1),'k-','LineWidth',5)

[~, i] = min(abs(x - Delta));
[~, j] = min(abs(x + Delta));
plot(x(i), U(i,k), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
plot(x(j), U(j,k), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
xlim([-pi,pi])
xlabel('x','FontSize',20);

nexttile
hold on
contour(B1, D1, F1, [0 0], 'k--', 'LineWidth', 5);
contour(B2, D2, F2, [0 0], 'k--', 'LineWidth', 5);
plot(betas2, Delta_sols2u, 'k','LineWidth',5);
plot(betas1, Delta_sols1u, 'k','LineWidth',5);

ylim([0,pi/2+0.1])
xlim([0,24])
set(gca,'XScale','log');
xlabel('\beta','FontSize',20);
ylabel('\Delta','FontSize',20);


% Creates the plots for Figure 7. Requires first running a finite difference simulation with a perturbation to the activity variable.

figure;
t = tiledlayout(2,3);

ax4 = nexttile;
plot(ax4, x, Q(:,70), '-','Color',[0 0.6 0.3 1], 'LineWidth',7.5); hold on;
plot(ax4,x, Q(:,1), ':','Color',[0 0.6 0.3 0.5],'LineWidth',7.5)
xlim([-pi,pi])
ylim([0.77,1.05])

ax5 = nexttile;
plot(ax5, x, Q(:,3000), '-','Color',[0 0.6 0.3 1], 'LineWidth',7.5); hold on;
plot(ax5,x, Q(:,1), ':','Color',[0 0.6 0.3 0.5],'LineWidth',7.5)
xlim([-pi,pi])
ylim([0.77,1.05])


ax6 = nexttile;
plot(ax6, x, Q(:,10000), '-','Color',[0 0.6 0.3 1], 'LineWidth',7.5); hold on;
plot(ax6,x, Q(:,1), ':','Color',[0 0.6 0.3 0.5],'LineWidth',7.5)
xlim([-pi,pi])
ylim([0.77,1.05])


ax7 = nexttile;
plot(ax7, x, A(:,70), '-','Color',[0.85 0.35 0.75 1], 'LineWidth',7.5); hold on;
plot(ax7,x, A(:,1), ':','Color',[0.85 0.35 0.75 0.5],'LineWidth',7.5)
xlim([-pi,pi])
ylim([0.092,0.108])


ax8 = nexttile;
plot(ax8, x, A(:,3000), '-','Color',[0.85 0.35 0.75 1], 'LineWidth',7.5); hold on;
plot(ax8,x, A(:,1), ':','Color',[0.85 0.35 0.75 0.5],'LineWidth',7.5)
xlim([-pi,pi])
ylim([0.092,0.108])


ax9 = nexttile;
plot(ax9, x, A(:,10000), '-','Color',[0.85 0.35 0.75 1], 'LineWidth',7.5); hold on;
plot(ax9,x, A(:,1), ':','Color',[0.85 0.35 0.75 0.5],'LineWidth',7.5)
xlim([-pi,pi])
ylim([0.092,0.108])

