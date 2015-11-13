r1 = linspace(0,3,100);
r2 = linspace(-3,0,100);
z = exp(4.*r);
y = exp(-4.*r);
plot(r1,y,'r')
hold('on')
plot(r2,z,'r')
xlabel('r','fontsize',15)
ylabel('exp(-4r)','fontsize',15)