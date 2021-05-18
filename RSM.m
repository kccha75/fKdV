% Function uses residual smoothing method (see Canuto)

% Inputs:
% v - initial guess
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient
% f - RHS
% numit - number of iterations

% Output:
% vnew - updated guess after iteration

% Other parameters
% alpha - weighting factors (see Canuto)
% beta - weighting factors (see Canuto)

% Note: pick either scaled weighted residual or weighted scaled residual
%--------------------------------------------------------------------------
function v=RSM(v,k,a,b,f,numit)
global L

if numit==0
    return;
end

alpha=0.380125;
beta=0.138155;

w=(L/length(v))^2;

for i=1:numit
    
    r=findR(f,v,k,a,b);

    % r(i) <- beta*r(i-1)+alpha*r(i)+beta*r(i+1) (see papers / Canuto)
    r=alpha*r+beta*[r(2:end);r(1)]+beta*[r(end);r(1:end-1)];
    
    v=v+w.*r;
    
end

end