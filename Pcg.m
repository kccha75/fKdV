% Function uses Preconditioned Conjugate Gradient to solve system Au=f

% Iterates until convergence condition met (small enough residual) OR
% maximum iterations reached

% Inputs:
% v - initial guess
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient
% f - RHS of equation
% H - preconditioner

% Ouputs:
% v - best estimate

% Optional display messages can be commented out or left in

function [v,r,i]=Pcg(v,k,a,b,f,H)

% Maximum number of iterations
numit=100;

r=f-LA(v,k,a,b);
d=H\r;
delta=r'*d;
if r==0
    return
end

for i=1:numit
    
    q=LA(d,k,a,b);
    alpha=delta/(d'*q);
    v=v+alpha*d; % estimate of new solution
    r=r-alpha*q;
%    disp(rms(r)) % display residual for reference
    
     if rms(r)<1e-12
%         fprintf('Preconditioned Conjugate Gradient Converged after %d iterations!\n',i);
        break
     end
    
    % update for next guess
    s=H\r;
    deltanew=r'*s;
    beta=deltanew/delta;
    d=s+beta*d;
    delta=deltanew;
    
end
if i==numit
    fprintf('Warning: Conjugate Gradient did not converge\n')
end
end