%{
n = input('n=');
a = zeros(n);

for i=1:n
    for j=1:n
        fprintf('a(%g,%g)=',i,j);
        a(i,j) = input('');
    end
end
%}
n = 5 ;
A = [1 2 3 5 6;
    4 25 6 7 2;
    7 48 -9 18 -5;
    10 11 12 15 -77;
    5 6 -4 -21 2];
AA = A;
H = A;
EPS = 1.e-16; % epsilon masina

% aducerea la forma superior Hessenberg

for k=1:n-2
    u = zeros(n,1);
    u(k+1) = A(k+1,k) + sign(A(k+1,k))*norm(A(k+1:n,k),2);
    u(k+2:n) = A(k+2:n,k);
    beta = u' * u / 2;
    if beta < EPS
        fprintf('beta nul sau foarte mic!\n');
        return
    end
    U = eye(n) - u * u' / beta
    % spre deosebire de agloritmul de la factorizarea QR, aici A se
    % inmulteste nu doar la stanga cu U, ci si la dreapta cu U transpus,
	% care este de fapt tot U, pentru ca U  este simetrica,
    % astfel pastrandu-se valorile proprii
    A = U * A * U;
end
% H este matrice ortogonal asemenea cu A
% ceea ce inseamna ca are aceleasi valori proprii ca A
% se poate verifica apeland eig(H) si eig(AA) !!
fprintf('Forma superior Hessenberg:\n');
H = A
input('');

%incep sa aplic pasi QR
flag = 1;
while (flag == 1)
    % etapa 1 - alegerea miu
    miu = H(n,n);
    % etapa 2 - se scade de pe diagonala princiapala deplasarea
    H = H - eye(n) * miu;
    % etapa 3 - triangularizarea ortogonala
    Q = eye(n);
    for k = 1: n-1
        % as putea folosi matrici Householder dar e vorba de a zeroriza un
        % singur element (celelalte sunt oricum deja zero); voi folosi, atunci,
        % matrici de rotatie plana Givens...
        rk = sqrt(H(k+1,k)^2 + H(k,k)^2);
        ck = H(k,k) / rk;
        dk = H(k+1,k) / rk;
        Rrond = [ck dk;-dk ck];
        Pk = eye(n);
        Pk(k:k+1,k:k+1) = Rrond;
        Q = Q * Pk';
        % voi inmulti pe H mai tarziu la dreapta cu el pt a pastra
        % valorile proprii
        H = Pk * H
    end;
    % la sfarsit, dupa transformari, matricea H devine R
    R = H;
    % etapa 4 - se construieste urmatorul element al sirului
    % H se inmulteste la stanga cu Qt si la dreapta cu Q, pentru a realiza
    % o TRANSFORMARE ORTOGONALA care PASTREAZA VALORILE PROPRII; 
    % H este Qt * H * Q, unde Qt * H este deja notat cu R aici
    H = R * Q;
    % etapa 5 - se adauga deplasarea pe diagonala principala
    H = H + miu * eye(n);
    % etapa 6 - testul de decuplare
    % se zerorizeaza fortat elementele de pe subdiagonala principala care
    % au scazut sub o precizie impusa
    for k=1:n-1
        if abs(H(k+1,k)) <= EPS * ( abs(H(k,k)) + abs(H(k+1,k+1)))
            H(k+1,k) = 0;
        end
    end
    H
    flag = 0;
    % etapa 7 - se aplica testele
    % testul 1
    for k=1:n-2
        if (H(k+1,k)~=0 && H(k+2,k+1)~=0)
            flag = 1;
        end
    end
    % testul 2
    if flag~=0
        k = 1;
        while (k<=n-1 && flag==0)
            a_ec = 1;
            b_ec = - H(k,k) - H(k+1,k+1);
            c_ec = H(k,k)*H(k+1,k+1) - H(k+1,k)*H(k,k+1);
            delta = b_ec^2 - 4*a_ec*c_ec;
            if delta>=0
                flag = 1;
            end
        end
    end
end
fprintf('Forma canonica Schur:\n');
H
input('');
% H este acum in forma canonica Schur;
% trebuie extrase valorile proprii
i = 1;
vp = zeros(1,n);
while i <= n-1
    if H(i+1,i) ~= 0
        % daca sub diagonala e nul, inseamna ca am vp reala fix pe
        % diagonala
        vp(i) = H(i,i);
        i = i + 1;
    else
        % altfel, inseamna ca am doua vp complex conjugate care sunt
        % vp ale blocului de ordin 2
        a_ec = 1;
        b_ec = - H(i,i) - H(i+1,i+1);
        c_ec = H(i,i)*H(i+1,i+1) - H(i+1,i)*H(i,i+1);
        delta = b_ec^2 - 4 * a_ec * c_ec;
        %vp(i:i+1) = roots([a_ec b_ec c_ec])
        vp(i) = (-b_ec + sqrt(-1) * sqrt(-delta) ) / (2 * a_ec);
        vp(i+1) = (-b_ec - sqrt(-1) * sqrt(-delta) ) / (2 * a_ec) ;
        i = i + 2;
    end
    if i==n
        vp(i) = H(n,n);
    end
end

vp