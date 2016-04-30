%APROXIMAREA sin CU SERII FOURIER

rel = 1;
while rel==1
    x=input('x=');
    n=input('Cati termeni? ');
    % calculez rezultatul folosind functia predefinita
    rez2 = sin(x);
    % aduc x in intervalul [0,2pi)
    while x<0
        x=x+2*pi;
    end
    while x>=2*pi
        x=x-2*pi;
    end
    minus = 0;                      % 1 daca in urma aducerii la primul cadran, rezultatul trebuie inmultit cu -1
    % aduc la primul cadran
    if (x>pi/2 && x<pi)             % cadranul II
        x = pi - x;
    elseif (x>=pi && x<3*pi/2)      % cadranul III
        x = x - pi;
        minus = 1;
    elseif (x>=3*pi/2 && x<2*pi)    % cadranul IV
        x = 2*pi-x; 
        minus = 1;
    end
    s = rem(n,2);
    % s = semnul alternant al termenilor dezvoltarii
    % 1 inseamna plus ; 0 inseamna minus
    rez = 0;
    for i= n:-1:1
        % pentru fiecare al i-lea termen calculez
        % x/1 * x/2 * ...  * x/(2n-1) 
        fact = 1;
        for j=2*i-1:-1:1
            fact = fact * (x/j);
        end
        % in functie de semn, adun sau scad termenul
        % dar evitand sa inmultesc termenul cu semnul,
        % adica sa calculez fact*s deoarece apar erori substantiale
        if s==1
           rez = rez + fact;
        else
           rez = rez - fact;
        end
        s=1-s;
    end
    if (s==1)
        rez = -rez;
    end
    if (minus==1)
        rez = -rez; 
    end
    fprintf('x=%g sinx=%g\n',rez,rez2);
    rel = input('Reluare? 1/0\n');
end
%end