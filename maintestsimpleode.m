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
L = 80;

% Number of grids (total points being 2^(grids))
grids=9;

% Grid where solution is exactly solved at 2^(gridexact) points
gridexact=8; 

% Iterations on down cycle (on finest grid and after restriction)
Nd=1; 

% Iterations on up cycle (after prolongation)
Nu=1; 

% Number of v-cycles
num_vcycles=10; 

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
b=-12*(sech(x+5).^2+sech(x-5).^2);
load('bnew.mat')
b=bnew;
% c(x) function
c=-3*ones(N,1);

A=LA(eye(N,N),k,a,b);
H=Hfd(L,N,a,b);

H2=Hfe(L,N,a,b);

RHS=sin(x);

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k, ...
    a,b,c,RHS,zeros(N,1));
cellH=setcellsH(vcyclegrid,H,L,cellN,cella,cellb);

% eig1=eig(cellH{1}\LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellb{1}));
% eig2=eig(cellH{2}\LA(eye(cellN{2},cellN{2}),cellk{2},cella{2},cellb{2}));
% eig3=eig(cellH{3}\LA(eye(cellN{3},cellN{3}),cellk{3},cella{3},cellb{3}));


[z,r]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellRHS,cellH,cellv);
disp(rms(r))

[z1,r1]=Pcg_it(zeros(N,1),k,a,b,RHS,H,20);
disp(rms(r1))

[z2,r2]=FMG_pre(pre_relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellRHS,cellH,cellv);
disp(rms(r2))

[z3,r3]=PMRR(zeros(N,1),k,a,b,RHS,H,200);
disp(rms(r3))