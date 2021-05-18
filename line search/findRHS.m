function RHS=findRHS(x,l,lambda)

% Define RHS function (2 bumps)
RHS=zeros(length(x),1);
for i=1:length(x)
    if abs(x(i)+l)<=1
        RHS(i)=cos(pi/2*(x(i)+l))^4;
    elseif abs(x(i)-l)<=1
        RHS(i)=cos(pi/2*(x(i)-l))^4;
    end
end
RHS=RHS*3*lambda;


end