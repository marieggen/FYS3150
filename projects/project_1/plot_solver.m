x = linspace(0,1,1000);
u = 1.0 - (1.0 - exp(-10.0)).*x - exp(-10.0.*x);

A = load('data_1000000.txt');
xx = A(:,1);
yy = A(:,2);

figure(1)
plot(xx , yy , 'r' , x , u , 'b')
legend('numerical','analytical');

doc = {'error_10.txt';'error_100.txt';'error_1000.txt';'error_10000.txt';'error_100000.txt';...
    'error_200000.txt';'error_300000.txt';'error_400000.txt';'error_500000.txt';...
    'error_600000.txt';'error_700000.txt';'error_800000.txt';'error_900000.txt';'error_1000000.txt'};
h = [];
relMaxError = [];

for i=1:10
    filename = doc{i};
    B = load(filename);
    relError = B(:,1);
    h = [ h B(1,2)];
    relMaxError = [relMaxError min(relError)];
end
figure(2)
plot(h,relMaxError)
