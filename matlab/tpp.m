% TRIANGULARIZAREA CU PIVOTARE PARTIALA - fara optimizari

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

for k=1:n
    % caut pivotul (nr cel mai mare IN MODUL pe coloana k)
    zona_cautare = abs(A(k:n,k));
    [el_max,poz_max] = max(zona_cautare);
    % corectie, pt ca poz_max e relativa la subcoloana in care caut
    poz_max = poz_max + k - 1;
    
    %{
    %VARIANTA 1, interschimbare manuala a liniilor
    if poz_max ~= k % daca am gasit maximul pe alta linie decat k
        aux = A(k,:); % interschimb ca sa aduc maximul pe linia k
        A(k,:) = A(poz_max,:);
        A(poz_max,:) = aux;
        aux = b(k);
        b(k) = b(poz_max);
        b(poz_max) = aux;
    end
    %}
    
    % VARIANTA 2 - se creeaza matricea P si se inmulteste
    P = eye(n);
    if poz_max~=k
        aux = P(k,:);
        P(k,:) = P(poz_max,:);
        P(poz_max,:) = aux;
        A = P * A;
        b = P * b;
    end
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

U = A
b

% substitutie inversa
x(n) = b(n) / A(n,n);
for i=n-1:-1:1
    sum = A(i,i+1:n) * x(i+1:n);
    x(i) = (b(i)-sum) / A(i,i);
end

x
    