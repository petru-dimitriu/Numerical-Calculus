% REZOLVARE A din R mxn IN SENSUL CELOR MAI MICI PATRATE

m = input('m=');
n = input('n=');
A = zeros(m,n);

for i=1:m
    fprintf('a(%g)=',i);
    A(i,:) = input('');
end

for i=1:m
    fprintf('b(%g)=',i);
    b(i) = input('');
end

Qtranspus = eye(m);
for k=1:n
    % vreau sa zerorizez coloana k
    % creez vectorul Householder u
    u = zeros(m,1);
    sigma = norm(A(k:m,k),2);
    u(k) = A(k,k) + sign(A(k,k))*sigma;
    u(k+1:m) = A(k+1:m,k);
    u
    % calculez beta
    beta = u' * u / 2
    % calculez matricea Householder
    U = eye(m) - u*u' / beta
    % inmultesc matricea A la stanga, ca sa zerorizez sub diagonala
    % principala
    A =  U * A
    % toate U-urile se acumuleaza intr-o matrice Qtranspus
    % Qt * A = R /, care inmultit la stanga cu Q da
    % Q * Qt * A = Q * R, iar Q * Qt = I_n pentru ca Q este ortogonala
    % A = QR, de unde vine si numele factorizare QR...
    Qtranspus = U * Qtranspus;
    % retin acest Qtranspus pt ca voi aplica aceleasi modificari vectorului
    % termenilor liberi
end

% la final, dupa toate transformarile, matricea A a devenit R
R = A
% aplicand aceleasi transformari si vectorului b, obtinem vectorul d
d = Qtranspus * b

%rezolv sistemul cu primele n ecuatii
x(n) = d(n,n) / R(n,n);
for i=n-1:-1:1
    sum = R(i,i+1:n)*d(i+1:n);
    x(i) = (b(i)-sum) / R(i,i);
end
% solutia in sensul cmmp...
x