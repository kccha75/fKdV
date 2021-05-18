% Function uses Preconditioned Minimum Residual Richardson iteration to
% solve system Au=f

% Inputs:
% v - best estimate
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient
% f - RHS of equation
% H - preconditioner
% numit - number of iterations

% Ouputs:
% v - best estimate after iteration

% See Boyd / Canuto for reference to algorithm

function [v,r]=PMRR(v,k,a,b,f,H,numit)

if numit==0
    return;
end

r=findR(f,v,k,a,b);
z=H\r;

if r==0
    return
end

for i=1:numit
    
    Az=LA(z,k,a,b);
    tau=dot(r,Az)/dot(Az,Az);
    v=v+tau*z;
    r=r-tau*Az;
    z=H\r;
    
end
end