% Function uses Conjugate Gradient to solve system Au=f as iterative
% smoother

% Iterates numit times as relaxation

% Inputs:
% v - initial guess
% k - wave number
% a - a(x) coefficient
% b - b(x) coefficient
% f - RHS of equation

% Ouputs:
% v - best estimate

% Optional display messages can be commented out or left in

function v=cg_it(v,k,a,b,f,numit)

if numit==0
    return;
end

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
    
%      if rms(r)<10e-15
% %         fprintf('Conjugate Gradient Converged after %d iterations!\n',i);
%         break
%      end
    
    % update for next guess
    deltanew=r'*r;
    beta=deltanew/delta;
    d=r+beta*d;
    delta=deltanew;
    
end

end