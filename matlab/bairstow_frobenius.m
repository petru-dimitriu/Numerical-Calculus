%1
reluareProgram=0;
while reluareProgram==0
    fprintf('introduceti ordinul polinomului ');
    n=input('n=');
 
    %2
    EPS = input('EPS=');
    a=zeros(1,n+1);
    aa=zeros(1,n+1);
    ac=zeros(1,n+1);
    %4
    reluare = 1;
    while reluare==1
        reluare = 0;
        for i=0:n
            ie=i+1;
            ip=n-i;
            fprintf('coeficient putere %g, a(%g)=', ip, ie);
            a(ie) = input('');
            aa(ie) = a(ie);
        end
        if a(1)==0
            fprintf('a(1) trebuie sa fie nenul\n');
            reluare = 1;
        end
    end
    aa1 = aa(1);
    for i=1:n+1
        aa(i) = aa(i)/aa1;
    end
 
    ac = a;
    nn =n;
    x = zeros(n,1);
    iter = 0;
    i = 0;
    b = zeros(n+3,1);
    c = zeros(n+2,1);
    while nn>2
        iter = iter + 1;
        fprintf ('Divizor nr %g\n',iter);
        b(1) = 0;
        b(2) = 0;
        c(1) = 0;
        c(2) = 0;
        p = 0; q = 0;
        flag = 1;
        while flag == 1
            for k =0:nn
                ibc = k+3;
                ia = k+1;
                b(ibc) = aa(ia) - p*b(ibc-1)- q*b(ibc-2);
                c(ibc) = b(ibc) - p*c(ibc-1)- q*c(ibc-2);
            end
            in = nn+3;
            di = c(in-2)^2-c(in-3)*(c(in-1)-b(in-1));
            pi = -b(in-1)*c(in-2)+b(in)*c(in-3);
            qi = -b(in) *c(in-2) + b(in-1)*(c(in-1) - b(in-1));
            if di == 0
                fprintf('di nul ');
                return;
            end
            p = p - pi/di;
            q = q - qi/di;
            rr = abs((pi+qi)/di);
            p
            q
            fprintf('r=\n%g',rr);
            if rr < EPS
                flag = 0;
            end
        end
        fprintf('Coeficienti divizor: \n');
        fprintf ('1, %g, %g\n',p,q);
        xd = roots([1 p q]);
        fprintf('Radacini divizor: %g\n',xd);
        i = i+1;
        x(i) = xd(1);
        i = i + 1;
        x(i) = xd(1);
        i = i + 1;
        x(i) = xd(2);
        nn = nn - 2;
        for k=0:nn
            ia = k+1;
            ib = k+3;
            aa(ia) = b(ib);
        end
    end
    iter = iter+1;
    fprintf('Divizor nr %g',iter);
    if nn==2
        fprintf('Coeficienti divizor:\n%g %g %g\n',aa(1),aa(2),aa(3));
        xd = roots([aa(1) aa(2) aa(3)]);
        i = i + 1;
        x(i) = xd(1);
        i = i + 1;
        x(i) = xd(2);
    else
        fprintf('Coeficienti divizor:\n%g %g\n',aa(1),aa(2));
        xd = roots([0 aa(1) aa(2)]);
        i = i + 1;
        x(i) = xd;
    end
    fprintf ('Radacini divizor:\n');
    xd
 
    % metoda Frobenius
    L1 = zeros(n-1,1);
    an = zeros(1,n);
    for i=0:n-1
        ir = n+1-i;
        an(i+1) = -a(ir)/a(1);
    end
    % bordare
    f = [L1 eye(n-1)];
    f = [f; an];
    %valori proprii matrice f
    xx = eig(f);
    %calcul radacini folosind roots
    x_matlab = roots(ac');
    %afisare comparativa
    fprintf('Radacini polinom - metoda Bairstow\n');
    x
    fprintf('Radacini polinom - metoda Frobenius\n');
    xx
    fprintf('Radacini polinom - functia roots\n');
    x_matlab
 
    fprintf('Doriti reluarea programului? \n 0=DA \n 1=NU \n');
    reluareProgram=input(' ');
end