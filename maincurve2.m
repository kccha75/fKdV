% clear;
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
tic
global L
% Domain size
L = 200;

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
b=ones(N,1);

% c(x) function
c=-3*ones(N,1);

% RHS function
l=10;

RHS=sech(x+l).^2+sech(x-l).^2;

% Preconditioner
preconditioner=@Hfd;

% Pseudo Arclength continuation step length
ds=0.01;

% Number of steps taken in pseudo arclength continuation
steps=20000;steps=230;

%--------------------------------------------------------------------------
% SOLVE FOR INITIAL SOLUTION + LAMBDA USING NEWTON
%--------------------------------------------------------------------------

% Initial lambda (guess for chosen gamma)
lambda0=0.1;

% Initial gamma (chosen)
gamma0=-8.2;

% Initial guess (near solution)
t=0.6;
v0=2*RHS;

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k, ...
    a,b,c,RHS,v0);

% Solve for initial solution using Newton SMG
[cellv,lambda0]=Newton_vcycle_tabletop_tol(lambda0,gamma0,x,pre_relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu, ...
    cella,cellb,cellc,cellRHS,preconditioner,cellv);

% Save solution
v0=cellv{1};

%--------------------------------------------------------------------------
% SOLVE FOR SECOND SOLUTION (NEAR INITIAL) USING NEWTON
%--------------------------------------------------------------------------
% Second lambda (guess for chosen gamma)
lambda1=0.15;
e=0.02;
% Second gamma (chosen)
gamma1=-8.3;

% % Initial guess uses previous solution (saved in cellv{1})
% A0=2*sech(x(1:N/2+1)+l).^2;
% A1=(1/6+1/6*1/2*(5*tanh(x(1:N/2+1)+l).^3-3*tanh(x(1:N/2+1)+l)));
% 
% b0=11/72;b2=115/504;b4=5/42;
% A2=-31/72/12+b0/12+b2/6*1/2*(3*tanh(x(1:N/2+1)+l).^2-1)-b4/8*1/8*(35*tanh(x(1:N/2+1)+l).^4-30*tanh(x(1:N/2+1)+l).^2+3);
% 
% y=A0+e*A1+e^2*A2;
% yy=[y;flipud(y(2:end-1))];
cellv{1}=v0;
% Solve for initial solution using Newton SMG
[cellv,lambda1]=Newton_vcycle_tabletop_tol(lambda1,gamma1,x,pre_relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu, ...
    cella,cellb,cellc,cellRHS,preconditioner,cellv);

% Save solution
v1=cellv{1};

%--------------------------------------------------------------------------
% FIND DV, DLAMBDA FOR PSEUDO ARCLENGTH INITIAL CONDITIONS
%--------------------------------------------------------------------------

lambda(2)=lambda1;
lambda(1)=lambda0;
gamma(2)=gamma1;
gamma(1)=gamma0;

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k, ...
    a,b,c,RHS,v0);

% Tangent approximation (found using secant approximation of two solutions)
dv=(v1-v0)/ds;
dlambda=(lambda(2)-lambda(1))/ds;
dgamma=(gamma(2)-gamma(1))/ds;

% Normalise
mag=sqrt(dot(dv,dv)+dlambda^2+dgamma^2);
dv=dv/mag;
dlambda=dlambda/mag;
dgamma=dgamma/mag;

% Pseudo arclength to trace solution family
[v,lambda,gamma]=pseudoarclength_test(v0,lambda0,gamma0,dv,dlambda,dgamma,ds,steps,x,pre_relaxation, ...
    vcyclegrid,cellN,cellk,Nd,Nu,cella,cellb,cellc,cellRHS,preconditioner);

toc