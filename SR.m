% Function uses Richardson Iteration to solve system Au=f

% Inputs:
% v - initial guess
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient
% f - RHS
% numit - number of iterations

% Output:
% vnew - updated guess after iteration

% ASSUMPTION min eigenvalue = 1, max eigenvalue = N^2/4
%--------------------------------------------------------------------------

function v=SR(v,k,a,b,f,numit)

w=2/(length(v)^2/4+1); % w - smoothing factor = 2/(eigmin+eigmax)
r=findR(f,v,k,a,b);
z=r;

for i=1:numit
    
    Az=LA(z,k,a,b);
    v=v+w*r;
    r=r-w*Az;
    z=r;
    
end

end