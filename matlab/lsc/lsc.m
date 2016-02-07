clear
clc

fprintf ('Input data in the form of a 2-column matrix containing data points\n')
input_data = input('input data = ');
format long e 

x = input_data(:,1);
y = input_data(:,2);
N = size(x,1);

xm = sum(x)/N; % mean value of x array
ym = sum(y)/N; 

u = x - xm; 
v = y - ym;

Suu = sum(u.*u);
Suv = sum(u.*v);
Suuu = sum(u.^3);
Suvv = sum(u.*v.*v);
Svv = sum(v.*v);
Svvv = sum(v.^3);
Svuu = sum(v.*u.*u);

A = [ Suu Suv ;
      Suv Svv ];

b = [1/2 * ( Suuu + Suvv) ;
    1/2 * (Svvv + Svuu) ];

sol = A\b;
uc = sol(1);
vc = sol(2);

% move found points, which were calculated in a different
% coordinate system, centered in (xm, ym)
xc = xm + uc;
yc = ym + vc;
% the centre of the circle is (xc, yc)
fprintf('Center of the LSC is (xc,yc)\n');
xc
yc

alfa = uc^2 + vc^2 + (Suu + Svv) / N;
R = sqrt(alfa);
fprintf('Radius of LSC:\n');
R

hold on
xlabel('Xlsc');
ylabel('Ylsc');
% plotting LSC
plot(xc,yc,'+r'); % plot the centre

for i=0:0.001:2*pi
    plot(xc+cos(i)*R,yc+sin(i)*R,'r','LineWidth',1.5);
end

% plot input data points
plot(x,y,'+b');