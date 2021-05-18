clear;clc
% 1D Fourier Newton Multigrid scheme on u_xx+au_x+bu+cu^2=f on domain 
% [-L/2,L/2] with periodic boundary conditions

% Uses vcycles to iteration (displays rms of residual after every v-cycle)

% Relaxation selections:
% SR - Stationary Richardson
% NSR - Non-Stationary Richardson (3 parameters)
% MRR - Minimum Residual Richardson (see Boyd)
% RSM - Residual smoothing method (see Canuto)
% cg - Conjugate Gradient

% Non-preconditioned relaxations ONLY!!!

% Important notes:
% SR and NSR assumes eigmin = 1, eigmax = N^2/4 (this is true for
% self-adjoint case only where a=1, b=0!!! with L=2*pi use with caution)

% Pseudo arclength continuation to trace solution family

%--------------------------------------------------------------------------
% INITIALIZE PARAMETERS
%--------------------------------------------------------------------------
tic
global L
% Domain size
L = 30;

% Number of grids (total points being 2^(grids))
grids=9;

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
relaxation=@MRR;

% Initial grid parameters
N = 2^grids; % initial grid size (Ideally 2^n sized)
k = 2*pi/L*[0:N/2-1 N/2 -N/2+1:-1]; % wave numbers
x = L*(-N/2:N/2-1)'/N; % Assumed to be periodic thus only going to N-1

% Initial lambda value
lambda0=1.5;
lambda1=1.4;

% a(x) function
a=zeros(N,1);

% b(x) function
b=6*ones(N,1);

% c(x) function
c=-9/2*ones(N,1);

% RHS function
l=L/10;
RHS=zeros(length(x),1);
for i=1:length(x)
    if abs(x(i)+l)<=1
        RHS(i)=cos(pi/2*(x(i)+l))^4;
    elseif abs(x(i)-l)<=1
        RHS(i)=cos(pi/2*(x(i)-l))^4;
    end
end
RHS=RHS*3;

% Pseudo Arclength continuation step length
ds=0.1;

% Number of steps taken in pseudo arclength continuation
steps=200;

%--------------------------------------------------------------------------
% SOLVE FOR INITIAL SOLUTION USING NEWTON
%--------------------------------------------------------------------------

% Initial guess (near solution)
v0=0.0*sech(x).^2;

% function to initiate cells
[cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k, ...
    a,b,c,RHS,v0);

% Select initial gamma value
lambda(1)=lambda0;

% Update variable wrt parameter
b1=lambda(1)*b;

% Set cells for new variable (required for parameters in a,b or c)
celltemp=setcellsNewton(vcyclegrid,cellb,b1);

% Solve for initial solution using Newton SMG
cellv=Newton_vcycle(relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu, ...
    cella,celltemp,cellc,cellRHS,cellv);

% Save solution
v1=cellv{1};

%--------------------------------------------------------------------------
% SOLVE FOR SECOND SOLUTION (NEAR INITIAL) USING NEWTON
%--------------------------------------------------------------------------

% Initial guess uses previous solution (saved in cellv{1})

% Select second gamma value (near initial)
lambda(2)=lambda1;

% Update variable for next gamma value
b2=lambda(2)*b;

% Set cells for new variable (required for parameters in a,b or c)
celltemp=setcellsNewton(vcyclegrid,cellb,b2);


% Solve for second solution using Newton SMG
cellv=Newton_vcycle(relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu, ...
    cella,celltemp,cellc,cellRHS,cellv);

% Save solution
v2=cellv{1};

%--------------------------------------------------------------------------
% FIND DV, DLAMBDA FOR PSEUDO ARCLENGTH INITIAL CONDITIONS
%--------------------------------------------------------------------------

% Tangent approximation (found using secant approximation of two solutions)
dv=(v2-v1)/ds;
dlambda=(lambda(2)-lambda(1))/ds;

% Normalise
mag=sqrt(dot(dv,dv)+dlambda^2);
dv=dv/mag;
dlambda=dlambda/mag;
     
% Pseudo arclength to trace solution family
[v,lambda]=pseudoarclength(v1,lambda,dv,dlambda,ds,steps,x,relaxation, ...
    vcyclegrid,cellN,cellk,Nd,Nu,cella,cellb,cellc,cellRHS);

toc

% Plot solution of A(-2) vs x
plot(lambda,v(round(-2*N/L+N/2),:),'x')
axis([0.5,3,-1,6])
xlabel('\lambda')
ylabel('A(-2)')