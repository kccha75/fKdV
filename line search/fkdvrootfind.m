clear;clc

% x points
L = 100;
grids=13;
N = 2^grids;
x = L*(-N/2:0)'/N;

% tolerance required
tol=1e-14;
odeopts = odeset('RelTol',1.e-14,'AbsTol',1.e-15);

% initial gamma values
gamma0=-1437.1;
gamma1=-1437.2;

ic=[0;0];
dc=0;

% solve for Ax at initial points
fkdv0 = @(x,y) [y(2);dc*y(1)-3*y(1)^2-gamma0*(sech(x).^2)];
[x0,y0]=ode45(fkdv0,x,ic,odeopts);
Ax0=y0(end,2);

fkdv1 = @(x,y) [y(2);dc*y(1)-3*y(1)^2-gamma1*(sech(x).^2)];
[x1,y1]=ode45(fkdv1,x,ic,odeopts);
Ax1=y1(end,2);

% loop secant method
while abs(Ax1)>tol
    
    gamma=(gamma0*Ax1-gamma1*Ax0)/(Ax1-Ax0);
    fkdv=@(x,y) [y(2);dc*y(1)-3*y(1)^2-gamma*(sech(x).^2)];
    [x,y]=ode45(fkdv,x,ic,odeopts);
    
    Ax=y(end,2);
    
    Ax0=Ax1;
    Ax1=Ax;
    
    gamma0=gamma1;
    gamma1=gamma;
    
end

x2=[x(1:end-1);flipud(-x(2:end))];
A2=[y(1:end-1,1);flipud(y(2:end,1))];

plot(x2,A2);xlabel('x');ylabel('A');title(['Solitary wave solution at \gamma=',num2str(gamma)])

gamma0=gamma;
v0=A2;