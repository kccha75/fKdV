% Function to find residual for 1D Fourier Multigrid
% r=f-A*v

% Inputs:
% f - RHS of Au=f
% v - best guess for solution
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient

% Outputs:
% r - residual for estimate v

function r=findR(f,v,k,a,b)

r=f-LA(v,k,a,b);

end