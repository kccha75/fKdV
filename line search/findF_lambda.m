% Function to find F_lambda (used in pseudo arclength)
% Partial derivative of F wrt lambda

% Inputs:
% x - vector of x ordinates
% gamma - parameter (not used)

% Ouputs:
% F_lambda - vector of partial derivative of F wrt lambda

function F_lambda=findF_lambda(a,b,c,u,x,l,lambda)

F_lambda=zeros(length(x),1);
for i=1:length(x)
    if abs(x(i)+l)<=1
        F_lambda(i)=-2*pi*cos(pi/2*(x(i)+l))^3*sin(pi/2*(x(i)+l));
    elseif abs(x(i)-l)<=1
        F_lambda(i)=2*pi*cos(pi/2*(x(i)-l))^3*sin(pi/2*(x(i)-l));
    end
end
F_lambda=-findRHS(x,l,1);

end