% TRIANGULARIZAREA SIMPLA - fara optimizari

%n = input('n = ');
n = 3;
A = zeros(n);
b = zeros(n,1);
x = zeros(n,1);
%{
for i=1:n
    for j=1:n
        fprintf('A(%g,%g)=',i,j);
        A(i,j) = input('');
    end
end

for i=1:n
    fprintf('B(%g)=',i);
    b(i) = input('');
end
%}

A = [4 6 3; 2 8 4; 7 3 1]; b = [2;43;5];
AA = A; bb = b;

for k=1:n
    % calculez ek
    ek = zeros(n,1);
    ek(k) = 1;
    % calculez vectorul de multiplicatori Gauss mk
    mk = zeros(n,1);
    for i=k+1:n
        mk(i) = A(i,k) / A(k,k);
    end
    % Mk = In - mk *ek transpus
    Mk = eye(n) - mk * ek'
    % am calculat Mk, il inmultesc la stanga lui A
    A = Mk * A
    % si la stanga lui b
    b = Mk * b;
end

% substitutie inversa
x(n) = b(n) / A(n,n);
for i=n-1:-1:1
    sum = A(i,i+1:n) * x(i+1:n);
    x(i) = (b(i)-sum) / A(i,i);
end

x
AA \ bb
