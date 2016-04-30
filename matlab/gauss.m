EPS = input('EPS=');
max_iter = input('max_iter=');
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
            fprintf('Rearanjati A');
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

nn = tril(a);
p = nn - a;
g = inv(nn)*p;
valp = eig(g);
ro = max(abs(valp));
ro
if ro<1
    fprintf('converge\n');
else
    fprintf('nu converge\n');
end


vninf = 1;
iter = 0;
format long e;
while (vninf>EPS) && (iter < max_iter)
    iter = iter + 1;
    xv = zeros(n,1);
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
        xn(i) = (b(i)-sum)/a(i,i);
    end
    vninf = max(abs(xn-xv));
    iter
    xn
    vninf
end
x = a \ b
format short;
iter
xn
x