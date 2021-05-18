clear;close all
% 1D Fourier Multigrid scheme on -u_xx+au_x+bu=f on domain [-L/2,L/2]
% with periodic boundary conditions

% Uses 1 FMG to iteration (displays rms of residual after every v-cycle)

% Relaxation selections:
% PSR - Stationary Richardson
% PNSR - Non-Stationary Richardson (3 parameters)
% PMRR - Minimum Residual Richardson (see Boyd)
% Pcg_it - Conjugate Gradient

% Preconditioned relaxations ONLY!!!

% Important notes:
% PSR and PNSR assumes eigmin = 1, eigmax = pi^2/4 with second-order FD
% preconditioning! (this is generally true)

%--------------------------------------------------------------------------
% INITIALIZE PARAMETERS
%--------------------------------------------------------------------------
global L

% Domain size
L = 20;

% Number of grids (total points being 2^(grids))
grids=10;

% Grid where solution is exactly solved at 2^(gridexact) points
gridexact=6; 

% Iterations on down cycle (on finest grid and after restriction)
Nd=1; 

% Iterations on up cycle (after prolongation)
Nu=1; 

% Number of v-cycles
num_vcycles=100; 

% Choice of relaxation
pre_relaxation=@Pcg_it;

% Initial grid parameters
N = 2^grids; % initial grid size (Ideally 2^n sized)
k = 2*pi/L*[0:N/2-1 N/2 -N/2+1:-1]; % wave numbers
x = L*(-N/2:N/2-1)'/N; % Assumed to be periodic thus only going to N-1

E1 = 1;
lambda = 1;
a1=2;

% RHS function
f=-1/2*a1*(-8-2*E1+lambda+(4+lambda)*cosh(2*x)).*sech(x).^4;

% a(x) function
a=zeros(N,1);

% b(x) function
b=-lambda+E1*(1-tanh(x).^2);

% Preconditioner
H=Hfd(L,N,a,b);

% Initial guess
% (not required since coarse grid guess is 0)

% Perform num_vcycles v-cycles
vcyclegrid=grids-gridexact+1; % number of vcycle grids 
% (since lower grids do not need storage as it is solved exactly)

% Initial guess (not needed for FMG)
v0=zeros(N,1);

% function to initiate cells
[cellN,cellk,cella,cellb,cellf,cellv]=setcells(vcyclegrid,N,k,a,b,f,v0);
cellH=setcellsH(vcyclegrid,H,L,cellN,cella,cellb);

%--------------------------------------------------------------------------
% SOLVE USING FMG + OPTIONAL VCYCLES AFTER
%--------------------------------------------------------------------------

% Perform Full Multigrid
% [cellv,r]=FMG_pre(pre_relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellb,cellf,cellH,cellv);
% fprintf('FMG\n %d\n',rms(r))

% V-cycles after
% [v,r]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellf,cellH,cellv);
% disp(rms(r))
% Speed testing

tic
for i=1:100
v1=Pcg_it(v0,k,a,b,f,H,8);
end
toc
% r1=findR(f,v1,k,a,b);
% fprintf('Preconditioned CG\n %d\n',rms(r1))
% 
% tic
% for i=1:100
% v2=LA(eye(N,N),k,a,b)\f;
% end
% toc
% r2=findR(f,v2,k,a,b);
% fprintf('Backslash\n %d\n',rms(r2))
% 
% tic
% for i=1:100
% v3=Pcg_it(v0,k,a,b,f,H,3);
% end
% toc
% r3=findR(f,v3,k,a,b);
% fprintf('PMRR\n %d\n',rms(r3))
% 
cellv{1}=v0;
tic
for i=1:100
[v4,r4]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellb,cellf,cellH,cellv);
end
toc
% fprintf('vcycle\n %d\n',rms(r4))