% GAUSS-SEIDEL + parametru de crestere a vitezei omega

EPS = input('EPS=');
nr_iteratii = input('nr max iteratii=');
pas = input('pasul variabilei omega=');
n = 0;
while n<2
    n = input('n=');
end
a = zeros(n);
ok = 0;
while (ok == 0)
    for i=1:n
        fprintf('Linia %g:',i);
        a(i,:)=input('');
    end
    ok = 1;
    for i=1:n
        if (abs(a(i,i))<EPS)
            fprintf('Elementul de pe diag nul sau foarte mic\n');
            ok = 0;
        end
    end
end

b = zeros(n,1);

for i=1:n
    fprintf('b(%g)=',i);
    b(i)=input('');
end

fprintf('Aproximatia initiala:\n');
xn = zeros(n,1);
for i=1:n
    fprintf('xn(%g) = ',i);
    xn(i)=input('');
end

nriter = zeros(1,fix(1/pas));
k = 0;

%ro_jacobi
nn = diag(diag(a));
p = nn - a;
g = inv(nn)*p;
valp = eig(g);
ro_jacobi = max(abs(valp));

%ro_gauss
nn = tril(a);
p = nn - a;
g = inv(nn)*p;
valp = eig(g);
ro = max(abs(valp));

%ro
if ro<1
    %fprintf('converge\n');
else
    fprintf('nu converge\n');
end
miniter = nr_iteratii+1;
for omega=1:pas:2 
    xn = zeros(n,1);
    k = k + 1;
    vninf = 1;
    iter = 0;
    format long e;
    xv = zeros(n,1);
    while (vninf>EPS) && (iter<nr_iteratii)
        iter = iter + 1;
        xv = xn;
        % adaptare gauss-seidel
        for i=1:n
            sum = 0;
            for j=1:i-1
                sum = sum + a(i,j)*xn(j);
            end
            for j=i+1:n
                sum = sum + a(i,j)*xv(j);   
            end
            xn(i) = (1-omega)*xv(i)+omega*((b(i)-sum)/a(i,i));
    %  xn(i) calculat cu metoda gauss --> |_________________|                 
        end
        vninf = max(abs(xn-xv));
        %iter
        %xn
        %vninf
    end
    %x = a \ b
    %format short;
    %iter
    %xn
    %x
    nriter(k)=iter;
	% pentru cazul cand numarul minim de iteratii se obtine pentru mai multe valori foarte apropiate ale omega
	% aleg drept omega_optim valoarea cea mai mare
    if (iter<=miniter)
		miniter=iter;
        omega_optim_aproximat_din_grafic=omega;
    end
        
        
end

format short;
omega_optim_aproximat_din_grafic
omega_optim_calculat = 2 / (1+sqrt(1-ro_jacobi^2))
ro_optim = ro_jacobi / (1+sqrt(1-ro_jacobi^2))
ro * omega_optim_calculat ^ 2/4

plot(1:pas:2,nriter);