clear all
close all

eigval_infile = ['eigval_2001.txt';'eigval_2002.txt';'eigval_2003.txt';'eigval_2004.txt'];
eigvec_infile = ['eigvec_2001.txt';'eigvec_2002.txt';'eigvec_2003.txt';'eigvec_2004.txt'];
omega_r = ['omega_r1.txt';'omega_r2.txt';'omega_r3.txt';'omega_r4.txt'];

for i=1:4
    A = load(eigval_infile(i,:));
    B = load(eigvec_infile(i,:));

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

%     save(omega_r(i),normEigvec(:,1));
    rho = linspace(1,10,n+2);
    plot(rho,normEigvec(:,1))
    hold('on')
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
