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
S_f = eye(n); % trebuie retinute permutarile de coloane, care se traduc in permutari in vectorul x

for k=1:n
    % caut pivotul (nr cel mai mare IN MODUL pe coloana k)
    zona_cautare = abs(A(k+1:n,k+1:n));
    [linie_max,poz_linie_max] = max(zona_cautare);
    [el_max,poz_coloana_max] = max(linie_max);
    % corectie, pt ca poz_max e relativa la subcoloana in care caut
    poz_linie_max = poz_linie_max + k - 1;
    poz_coloana_max = poz_coloana_max + k - 1;
    
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
    S = eye(n);
    if poz_linie_max~=k
        aux = P(k,:);
        P(k,:) = P(poz_linie_max,:);
        P(poz_linie_max,:) = aux;
        A = P * A
        b = P * b
    end
    if poz_coloana_max~=k
        aux = S(:,k);
        S(:,k) = S(:,poz_coloana_max);
        S(:,poz_coloana_max) = aux;
        A = A * S;
        S_f = S_f * S
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

% substitutie inversa
x(n) = b(n) / A(n,n);
for i=n-1:-1:1
    sum = A(i,i+1:n) * x(i+1:n);
    x(i) = (b(i)-sum) / A(i,i);
end

x = S_f * x;

x
    