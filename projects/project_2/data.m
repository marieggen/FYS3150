clear all
close all
format long

% eigval_infile = ['eigval_2001.txt';'eigval_2002.txt';'eigval_2003.txt';'eigval_2004.txt'];
% eigvec_infile = ['eigvec_2001.txt';'eigvec_2002.txt';'eigvec_2003.txt';'eigvec_2004.txt'];

% eigval_infile = 'eigval_300100.txt';
% eigvec_infile = 'eigvec_300100.txt';

% eigval_infile = ['eigval_2003.txt';'eigval_2004.txt'];
% eigvec_infile = ['eigvec_2003.txt';'eigvec_2004.txt'];

% eigval_infile ='eigval_2002.txt';
% eigvec_infile ='eigvec_2002.txt';

eigval_infile ='eigval_2001.txt';
eigvec_infile ='eigvec_2001.txt';

rho_file = 'rho.txt';


%for i=1:4
for i =1:1
    A = load(eigval_infile(i,:));
    B = load(eigvec_infile(i,:));
    rho = load(rho_file(:));

    n = B(1);
    eigval_rand = A(:,1);
    eigvecT_rand = [];
    m = n+1;
    k = 2;

    for i=1:n
        eigvecT_rand = [eigvecT_rand, B(k:m)];
        k = m+1;
        m = m + n;
    end
    eigvec_rand = eigvecT_rand';

    [eigval,index] = sort(eigval_rand);

    eigvecTemp = [];
    for i=1:n
        eigvecTemp = [eigvecTemp, eigvec_rand(:,index(i))];
    end

    normEigvec = [];
    for i=1:n
        normalize = eigvecTemp(:,i)./norm(eigvecTemp(:,i));
        boundaryEigvec = [0 ; normalize ; 0];
        normEigvec = [normEigvec, boundaryEigvec];
    end
    
    prob = normEigvec.^2;
    rho_min = rho(1,1);
    rho_max = rho(1,2);
    
    %normalization of probability
   rho = linspace(rho_min,rho_max,n+2);
   h = rho(2)-rho(1);
   vec_1 = prob(:,1);
   %describe analytically
   C = 1.0/(sum(vec_1)*h);
   vec_1C = C*vec_1;
    
    plot(rho,vec_1C)
    xlabel('\rho = r/\alpha [dim.less]','fontsize',16)
    ylabel('|\psi(\rho)|^2 [dim.less]','fontsize',16)
    h = legend('\omega_r = 0.01','\omega_r = 5.0','fontsize',20);
    set(h,'FontSize',16)
    hold('on')
    
    %print three first eigenvalues
%     eigval(1)
%     eigval(2)
%     eigval(3)
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
