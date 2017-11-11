%=========================================================================
% I.Y. Tanasa (1034117)       
% Aleman Zapata, R.A. (59383)
% Modeling Dynamics Project
%=========================================================================
%% QUESTION 7 VERIVY THE POSITIVE INVARIANT FOR UNDRIVEN SYSTEM
% close all, clear all, clc;
syms e f p q u x y z
% e= 10^-2; f= 0.2; p= 0.5; q= 0.05; % condition 1
e= 10^-4; f= 1; p= 100; q= 10^-6; % condition 2
u=0 ; % without external excitation
eq11= (x + y - q*x^2 - x*y + u) / e == 0;% 1st order system x
eq12= -y + f*z - x*y + u == 0; % 1st order system y
eq13= (x-z) / p == 0; % 1st order system z
solxyz2= solve([eq11, eq12, eq13 ], [x,y,z]); % 3 set of fp (x,y,z)
fp1q7= [double(solxyz2.x(1,1)) double(solxyz2.y(1,1))...
    double(solxyz2.z(1,1))];
fp2q7= [double(solxyz2.x(2,1)) double(solxyz2.y(2,1))...
    double(solxyz2.z(2,1))];
fp3q7= [double(solxyz2.x(3,1)) double(solxyz2.y(3,1))...
    double(solxyz2.z(3,1))];

return
%%
func = @(t,x) [(x(1) + x(2) - q*x(1).^2 - x(1)*x(2) ) / e; -x(2) + f*x(3) - x(1)*x(2) ; (x(1)-x(3)) / p];
% [t,xa] = ode45(func,[0 100],[0.5 0.5 0.5]);
[t,xa] = ode45(func,[0 .10],[1400    1.5    1400]);

plot3(xa(:,1),xa(:,2),xa(:,3));
hold on

var1 =[ 0 0 0];
var2 =[ 1413.71          1.00       1413.71];

% 
% var2 =[ 19/2 - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f,...
%     (3*f)/4 + (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80, 19/2 ...
%     - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f];
% var3 =[ (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2, (3*f)/4 ...
%     - (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80,...
%     (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2];


plot3(var1(1),var1(2),var1(3),'r*')
plot3(var2(1),var2(2),var2(3),'g*')

% plot3(var1(1),var1(2),var1(3),'r*')
% plot3(var2(1),var2(2),var2(3),'g*')

%plot3(var3(1),var3(2),var3(3),'b*')

grid on
title('Solution curve')
%%


% attempt on quiver
x= -3:0.5:3;
y= -3:0.5:3;
z= -3:0.5:3;
[X,Y,Z] = meshgrid(x, y, z);
e= 10^-4; f= 1; p= 100; q= 10^-6; % condition 2
u=0 ; % without external excitation
eq11= (x + y - q.*x.^2 - x.*y + u) / e == 0;% 1st order system x
eq12= -y + f.*z - x.*y + u == 0; % 1st order system y
eq13= (x-z) / p == 0; % 1st order system z
% using simplified equation
h= 100.*x - 5.*x.^2 + (100.*f.*x)/(x + 1) - (100.*f.*x.^2)/(x + 1);
i= 100.*x - 5.*x.^2 + (100.*f.*x)/(x + 1) - (100.*f.*x.^2)/(x + 1);
j= 100.*x - 5.*x.^2 + (100.*f.*x)/(x + 1) - (100.*f.*x.^2)/(x + 1);
[R,S,T] = meshgrid(eq11, eq12, eq13);
[H,I,J] = meshgrid(h, i, j);
figure
quiver3(X,Y,Z,R,S,T)
view(-35,45)
figure
quiver3(X,Y,Z,H,I,J)
view(-35,45)

