clear;close
% 1D Fourier Newton Multigrid scheme on u_xx+au_x+bu+cu^2=f on domain 
% [-L/2,L/2] with periodic boundary conditions

% Uses vcycles to iteration (displays rms of residual after every v-cycle)

% Relaxation selections:
% SR - Stationary Richardson
% NSR - Non-Stationary Richardson (3 parameters)
% MRR - Minimum Residual Richardson (see Boyd)
% RSM - Residual smoothing method (see Canuto)
% cg_it - Conjugate Gradient

% Non-preconditioned relaxations ONLY!!!

% Important notes:
% SR and NSR assumes eigmin = 1, eigmax = N^2/4 (this is true for
% self-adjoint case only where a=1, b=0!!! with L=2*pi use with caution)

%--------------------------------------------------------------------------
global L

% Domain size
L = 30;

% Number of grids (total points being 2^(grids))
grids=9;

% Grid where solution is exactly solved at 2^(gridexact) points
gridexact=5; 

% Iterations on down cycle (on finest grid and after restriction)
Nd=3; 

% Iterations on up cycle (after prolongation)
Nu=1; 

% Number of v-cycles
num_vcycles=1; 

% Choice of relaxation
relaxation=@MRR;

% Initial grid parameters
N = 2^grids; % initial grid size (Ideally 2^n sized)
k = 2*pi/L*[0:N/2-1 N/2 -N/2+1:-1]; % wave numbers
x = L*(-N/2:N/2-1)'/N; % Assumed to be periodic thus only going to N-1

% Define RHS as in thesis
% Equation Parameters (see notes)
l=0.5;
RHS=zeros(N,1);
for i=1:length(x)
    if abs(x(i))<=l
        RHS(i)=(4/(3*l))*cos(pi/2*x(i)./l)^4;
    else
        RHS(i)=0;
    end
end
gamma=-1;
RHS=RHS*gamma;

% a(x) function
a=zeros(N,1);

% b(x) function
b=1.7454*ones(N,1);

% c(x) function
c=-9/2*ones(N,1);

% Initial guess (near solution)
v0=0.0*sech(x).^2;

% Perform num_vcycles v-cycles
vcyclegrid=grids-gridexact+1; % number of vcycle grids 
% (since lower grids do not need storage as it is solved exactly)

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k,a,b,c,RHS,v0);

% tic
cellv=Newton_vcycle(relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellc,cellRHS,cellv);
% toc

% Plot solution
% fig1=figure;
% plot(x,v);title('Solution');