x = linspace(0,1,1000);
u = 1.0 - (1.0 - exp(-10.0)).*x - exp(-10.0.*x);

A = load('data_1000.txt');
xx = A(:,1);
yy = A(:,2);

figure(1)
plot(xx , yy , 'r' , x , u , 'b')
legend('numerical','analytical');
title('Comparison; numerical and analytical solution','FontSize',20);
xlabel('x','FontSize',20);
ylabel('u(x)','FontSize',20)
set(gca,'FontSize',15);

doc = {'error_10.txt';'error_100.txt';'error_1000.txt';'error_10000.txt';'error_100000.txt';...
    'error_200000.txt';'error_300000.txt';'error_400000.txt';'error_500000.txt';...
    'error_600000.txt';'error_700000.txt';'error_800000.txt';'error_900000.txt';'error_1000000.txt'};

B = load('time_lu_1000.txt');
y = B(:,1);
y = [0;y;0];
figure(3)
plot(xx,y,'r')
title('Solution with LU-decomposition','FontSize',20);
xlabel('x','FontSize',20);
ylabel('u(x)','FontSize',20)
set(gca,'FontSize',15);
% doc = {'error_10.txt';'error_100.txt'};

h = [];
relMaxError = [];

for i=1:14
    filename = doc{i};
    B = load(filename);
    relError = B(:,1);
    h = [ h B(1,2)];
    relMaxError = [relMaxError min(relError)];
end
figure(2)
plot(h,relMaxError,'r')
title('Error analysis','FontSize',20);
ylabel('log10(Maximum relative error)','FontSize',20);
xlabel('log10(step length)','FontSize',20);
set(gca,'FontSize',15);

















