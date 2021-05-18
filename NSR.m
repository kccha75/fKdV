% Function uses NSR to solve system Au=f
% Non stationary Richardson 3 parameters

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

% ASSUMPTION min eigenvalue = 1, max eigenvalue = N^2/4
%--------------------------------------------------------------------------

function v=NSR(v,k,a,b,f,numit)

if numit==0
    return;
end

l=3;
w=zeros(l,1);
for j=1:l
    w(j)=((length(v)^2+4)/8+(length(v)^2-4)/8*cos((2*j-1)*pi/(2*l)))^(-1); % w - smoothing factor for preconditioned NSR (see Canuto)
end

r=findR(f,v,k,a,b);
z=r;

for i=1:numit
    
    Az=LA(z,k,a,b);
    v=v+w(mod(i,3)+1)*z;
    r=r-w(mod(i,3)+1)*Az;
    z=r;
    
end

end