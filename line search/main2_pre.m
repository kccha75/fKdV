clear;close
% 1D Fourier Newton Multigrid scheme on u_xx+au_x+bu+cu^2=f on domain 
% [-L/2,L/2] with periodic boundary conditions

% Uses vcycles to iteration (displays rms of residual after every v-cycle)

% Relaxation selections:
% PSR - Stationary Richardson
% PNSR - Non-Stationary Richardson (3 parameters)
% PMRR - Minimum Residual Richardson (see Boyd)
% PRSM - Residual smoothing method (see Canuto)
% Pcg_it - Conjugate Gradient

% Preconditioned relaxations ONLY!!!

% Important notes:
% PSR and PNSR assumes eigmin = 1, eigmax = pi^2/4 with second-order FD
% preconditioning! (this is generally true)

%--------------------------------------------------------------------------
global L

% Domain size
L = 100;

% Number of grids (total points being 2^(grids))
grids=13;

% Grid where solution is exactly solved at 2^(gridexact) points
gridexact=5; 

% Iterations on down cycle (on finest grid and after restriction)
Nd=1; 

% Iterations on up cycle (after prolongation)
Nu=1; 

% Number of v-cycles
num_vcycles=1; 

% Choice of relaxation
pre_relaxation=@PMRR;

% Initial grid parameters
N = 2^grids; % initial grid size (Ideally 2^n sized)
k = 2*pi/L*[0:N/2-1 N/2 -N/2+1:-1]; % wave numbers
x = L*(-N/2:N/2-1)'/N; % Assumed to be periodic thus only going to N-1

% RHS function
l=10;
gamma=-8;
RHS=(sech(x+l).^2+sech(x-l).^2);%RHS=sech(x).^2;

% a(x) function
a=zeros(N,1);

% b(x) function
b=0*ones(N,1);

% c(x) function
c=-3*ones(N,1);

% Preconditioner
preconditioner=@Hfd;

% Initial guess (near solution)
v0=1.9*RHS;

% v0=RHS.*[zeros(N/2,1);ones(N/2,1)];
% Perform num_vcycles v-cycles
vcyclegrid=grids-gridexact+1; % number of vcycle grids 
% (since lower grids do not need storage as it is solved exactly)

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k,a,b,c,RHS,v0);

% tic
cellv=Newton_vcycle_pre(gamma,pre_relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellc,cellRHS,preconditioner,cellv);
% toc