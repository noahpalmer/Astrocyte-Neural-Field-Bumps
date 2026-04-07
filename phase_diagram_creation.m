% This file can be used to recreate the main numerical results of the paper (i.e. Figure 5 - 7). 
% The first section contains code needed to create Figures 5 and 6 and the second section for Figure 7.

tol = 1e-11; % For determing whether small eigenvalues in "fourierphase.m" are positive or zero.
bumpwidth = 1.5; % Initial guess for bump width, setup for larger bump.


% Values in parameter space for analytic phase diagrams of "fourierphase.m" and "analyticphase.m". Numerical simulations are quite slow for high resolutions
% and so code is currently setup for lower resolution there.
beta_vals = linspace(0.001,0.1,300);
D_vals = linspace(0.001,0.6,300);
gamma_vals = linspace(1,1.6,300);

% Values in parameter space for "numericalphase.m". Setup to be lower resolution than the analytic comps for speed.
beta_valsn = linspace(0.001,0.1,50);
D_valsn = linspace(0.001,0.6,50);
gamma_valsn = linspace(1,1.6,50);

% Parameters for numerical sims
N = 8000; % Number of spatial grid points
T = 500; % Simulation length, function has built in time step length of 0.01.

% Other parameters
D = 1;
gamma = 2;
theta = 0.1;

% Computes a phase diagram using the Fourier truncation method. Returns a
% binary matrix with 1s for unstable locations in parameters space and 0s
% for stable locations. Last argument of the function tells it whether the
% phase diagram is (D,beta): 1 or (gamma,beta): 2.

stablevalsD = fourierphase(beta_vals,D_vals,gamma,theta,tol);
stablevalsgamma = fourierphase(beta_vals,D,gamma_vals,theta,tol);


% Computes the boundary of the unstable region by approximating the
% derivative of the Evans function at lambda=0 and determining where it is 
% 0. Last argument of the function tells it whether we are considering
% finite D: 1 or if we're considering the large diffusion limit: 2.

phaseD = analyticphase(beta_vals,D_vals,gamma,theta,1);
phasegamma = analyticphase(beta_vals,D,gamma_vals,theta,1);

% Runs a simulation of the system with initial condition of the stationary
% bump solution, perturbed by the right shift perturbation in Eq. 5.2 for
% the given parameters. Computes the velocity of the perturbed bump at each location in phase space 
% at time T as well as the cumulative drift at time T. 

% (D,\beta) space diagrams

[velocityD, driftD] = numericalphase(beta_valsn,D_valsn,gamma,theta,N,T);

% (\gamma,\beta) space diagrams

[velocityg, driftg] = numericalphase(beta_valsn,D,gamma_valsn,theta,N,T);


% Figures

figure;
imagesc(D_valsn,beta_valsn,driftD); hold on
set(gca,'YDir','normal')
colormap([linspace(1,0,256)' linspace(1,0,256)' ones(256,1)])
colorbar;
contour(D_vals,beta_vals,phaseD,[0 0],'k-','LineWidth',3);
contour(D_vals,beta_vals,stabelvalsD,[0.1 0.1],'k--','LineWidth',3);
ylim([0.001,0.1])

figure;
imagesc(D_valsn,beta_valsn,velocityD); hold on
set(gca,'YDir','normal')
colormap([ones(256,1) linspace(1,0,256)' linspace(1,0,256)'])  
colorbar;
contour(D_vals,beta_vals,phaseD,[0 0],'k-','LineWidth',3);
contour(D_vals,beta_vals,stabelvalsD,[0.1 0.1],'k--','LineWidth',3);
ylim([0.001,0.1])

figure;
imagesc(gamma_valsn,beta_valsn,driftg); hold on
set(gca,'YDir','normal')
colormap([linspace(1,0,256)' linspace(1,0,256)' ones(256,1)])
colorbar;
contour(gamma_vals,beta_vals,phasegamma,[0 0],'k-','LineWidth',3);
contour(gamma_vals,beta_vals,stablevalsgamma,[0.1 0.1],'k--','LineWidth',3);
ylim([0.001,0.1])

figure;
imagesc(gamma_valsn,beta_valsn,velocityg); hold on
set(gca,'YDir','normal')
colormap([ones(256,1) linspace(1,0,256)' linspace(1,0,256)'])  
colorbar;
contour(gamma_vals,beta_vals,phasegamma,[0 0],'k-','LineWidth',3);
contour(gamma_vals,beta_vals,stabelvalsgamma,[0.1 0.1],'k--','LineWidth',3);
ylim([0.001,0.1])



%%

% Contains code for the D\to\infty limit phase diagram. 

% For the analytic boundaries of the unstable region for various D values
phaseD1 = analyticphase(beta_vals,1,gamma_vals,theta,1);
phaseD2andahalf = analyticphase(beta_vals,2.5,gamma_vals,theta,1);
phaseD10 = analyticphase(beta_vals,10,gamma_vals,theta,1);
phaseD50 = analyticphase(beta_vals,50,gamma_vals,theta,1);

% For the analytic boundary of the unstable region for D\to\infty limit
phaseDinf = analyticphase(beta_vals,1,gamma_vals,theta,2);  % D value does not matter if choice = 2

% Approximation of the velocity of stationary bumps at T for the D\to\infty limit

[velocityDinf, ~] = numericalphase(beta_valsn,1000,gamma_valsn,theta,N,T); % D is set to 1000 to simulate large diffusion.


figure;
imagesc(gamma_valsn,beta_valsn,velocityDinf); hold on
set(gca,'YDir','normal')
colormap([ones(256,1) linspace(1,0,256)' linspace(1,0,256)'])  
colorbar;
contour(gamma_vals,beta_vals,phaseD1,[0 0],'k--','LineWidth',3);
contour(gamma_vals,beta_vals,phaseD2andahalf,[0 0],'k--','LineWidth',3);
contour(gamma_vals,beta_vals,phaseD10,[0 0],'k--','LineWidth',3);
contour(gamma_vals,beta_vals,phaseD50,[0 0],'k--','LineWidth',3);
contour(gamma_vals,beta_vals,phaseDinf,[0 0],'k-','LineWidth',3);
ylim([0.001,0.1])