clear;
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
% preconditioning! (this is generally true for a~1)

% Pseudo arclength continuation to trace solution family

%--------------------------------------------------------------------------
% INITIALIZE PARAMETERS
%--------------------------------------------------------------------------

global L
% Domain size
L = 100;

% Number of grids (total points being 2^(grids))
grids=13;

% Grid where solution is exactly solved at 2^(gridexact) points
gridexact=8; 

% Iterations on down cycle (on finest grid and after restriction)
Nd=1; 

% Iterations on up cycle (after prolongation)
Nu=1; 

% Number of v-cycles
num_vcycles=1; 

% Perform num_vcycles v-cycles
vcyclegrid=grids-gridexact+1; % number of vcycle grids 
% (since lower grids do not need storage as it is solved exactly)

% Choice of relaxation
pre_relaxation=@PMRR;

% Initial grid parameters
N = 2^grids; % initial grid size (Ideally 2^n sized)
k = 2*pi/L*[0:N/2-1 -N/2 -N/2+1:-1]; % wave numbers
x = L*(-N/2:N/2-1)'/N; % Assumed to be periodic thus only going to N-1

% a(x) function
a=zeros(N,1);

% b(x) function
b=0*ones(N,1); % fixed parameter

% c(x) function
c=-3*ones(N,1);

% RHS function
l=10;

RHS=sech(x+l).^2+sech(x-l).^2;RHS=sech(x).^2;

% Preconditioner
preconditioner=@Hfd;

% Pseudo Arclength continuation step length
ds=0.05;

% Number of steps taken in pseudo arclength continuation
steps=5000;

%--------------------------------------------------------------------------
% SOLVE FOR INITIAL SOLUTION USING NEWTON
%--------------------------------------------------------------------------

% Initial gamma 
gamma0=0;

% Initial guess (near solution)
v0=0*RHS;

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k, ...
    a,b,c,RHS,v0);


% Solve for initial solution using Newton SMG
cellv=Newton_vcycle_pre(gamma0,pre_relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu, ...
    cella,cellb,cellc,cellRHS,preconditioner,cellv);

% Save solution
v1=cellv{1};

%--------------------------------------------------------------------------
% SOLVE FOR SECOND SOLUTION (NEAR INITIAL) USING NEWTON
%--------------------------------------------------------------------------

% Second gamma
gamma1=-.001;

% Initial guess uses previous solution (saved in cellv{1})

% Solve for second solution using Newton SMG
cellv=Newton_vcycle_pre(gamma1,pre_relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu, ...
    cella,cellb,cellc,cellRHS,preconditioner,cellv);

% Save solution
v2=cellv{1};

%--------------------------------------------------------------------------
% FIND DV, DLAMBDA FOR PSEUDO ARCLENGTH INITIAL CONDITIONS
%--------------------------------------------------------------------------

gamma(2)=gamma1;
gamma(1)=gamma0;

% Tangent approximation (found using secant approximation of two solutions)
dv=(v2-v1)/ds;
dgamma=(gamma(2)-gamma(1))/ds;

% Normalise
mag=sqrt(dot(dv,dv)+dgamma^2);
dv=dv/mag;
dgamma=dgamma/mag;
     
% Pseudo arclength to trace solution family
[v,gamma]=pseudoarclength_pre(v1,gamma,dv,dgamma,ds,steps,pre_relaxation, ...
    vcyclegrid,cellN,cellk,Nd,Nu,cella,cellb,cellc,cellRHS,preconditioner);

