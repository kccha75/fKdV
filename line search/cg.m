% Function uses Conjugate Gradient to solve system Au=f

% Iterates until convergence condition met (small enough residual) OR
% maximum iterations reached

% Inputs:
% v - initial guess
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient
% f - RHS of equation

% Ouputs:
% v - best estimate

% Optional display messages can be commented out or left in

function v=cg(v,k,a,b,f)

% Maximum number of iterations
numit=1000;

r=f-LA(v,k,a,b);
d=r;
delta=r'*r;
if r==0
    return
end

for i=1:numit
    
    q=LA(d,k,a,b);
    alpha=delta/(d'*q);
    v=v+alpha*d; % estimate of new solution
    r=r-alpha*q;
%    disp(rms(r)) % display residual for reference
    
     if rms(r)<10e-13
         fprintf('Conjugate Gradient Converged after %d iterations!\n',i);
        break
     end
    
    % update for next guess
    deltanew=r'*r;
    beta=deltanew/delta;
    d=r+beta*d;
    delta=deltanew;
    
end
if i==numit
    fprintf('Warning: Conjugate Gradient did not converge after %d iterations!\n',numit)
end