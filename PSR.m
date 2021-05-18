% Function uses Preconditioned Richardson Iteration to solve system Au=f

% Inputs:
% v - initial guess
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient
% f - RHS
% H - preconditioner
% numit - number of iterations

% Output:
% vnew - updated guess after iteration

% ASSUMPTION min eigenvalue = 1, max eigenvalue = pi^2/4 AFTER FD
% PRECONDITIONING
%--------------------------------------------------------------------------

function v=PSR(v,k,a,b,f,H,numit)

if numit==0
    return;
end

w=8/(4+pi^2); % w - smoothing factor for preconditioned SR (see Boyd)
r=findR(f,v,k,a,b);
z=H\r;

for i=1:numit
    
    Az=LA(z,k,a,b);
    v=v+w*z;
    r=r-w*Az;
    z=H\r;
    
end

end