plot(gamma,lambda)
xlabel('gamma');ylabel('lambda');title('Gamma vs Lambda for tabletop solutions')
figure;
plot(x,v(:,1:10:end))
title('Tabletop solution u vs x')
xlabel('x');ylabel('u')