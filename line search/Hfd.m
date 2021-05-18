% Function generates preconditioning matrix of operator A using second
% order finite differencing

% second order FD on -u_xx+au_x+bu

% Inputs:
% L - domain size
% N - grid size
% a - a(x) array
% b - b(x) array

% Outputs:
% H - preconditioning matrix (sparse)

function H=Hfd(L,N,a,b)

% Step size
dx=L/N;

% Preset diagonal array
B=zeros(N,3);

B(1:N,2)=(2+b(1:N)*dx^2)/dx^2;
B(1:N-1,1)=(-2-a(2:N)*dx)/(2*dx^2);
B(2:N,3)=(-2+a(1:N-1)*dx)/(2*dx^2);

% Diagonal position array
d=[-1 0 1];

% Generate sparse matrix
H=spdiags(B,d,N,N);

% First row (cyclic)
H(1,N)=(-2-a(1)*dx)/(2*dx^2);

% Last row (cyclic)
H(N,1)=(-2+a(N)*dx)/(2*dx^2);

end