xf = 1000;
gamma = linspace(-500,0,5000);
odeopts = odeset('RelTol',1.e-6,'AbsTol',1.e-12);
for ig = 1:length(gamma)
    dc = .0;
    gc = gamma(ig);
    fkdv = @(x,y) [y(2);dc*y(1)-3*y(1)^2-gc*(1-tanh(x).^2)];
    ic = 0.5*[-2/xf^2;4/xf^3];
    [xo,yo] = ode45(fkdv,[xf .0],ic,odeopts);
    af(ig) = yo(end,2);
end
clf
hold on
plot(gamma,af,'b')
plot(gamma,0*gamma,'-.')
