%=========================================================================
% I.Y. Tanasa (1034117)       
% Aleman Zapata, R.A. (59383)
% Modeling Dynamics Project for question 7 to 10
%=========================================================================
%% QUESTION 7 VERIVY THE POSITIVE INVARIANT FOR UNDRIVEN SYSTEM
% close all, clear all, clc;
% syms e f p q u x y z
% % e= 10^-2; f= 0.2; p= 0.5; q= 0.05; % condition 1
% e= 10^-4; f= 1; p= 100; q= 10^-6; % condition 2
% u=0 ; % without external excitation
% eq11= (x + y - q*x^2 - x*y + u) / e == 0;% 1st order system x
% eq12= -y + f*z - x*y + u == 0; % 1st order system y
% eq13= (x-z) / p == 0; % 1st order system z
% solxyz2= solve([eq11, eq12, eq13 ], [x,y,z]); % 3 set of fp (x,y,z)
% fp1q7= [double(solxyz2.x(1,1)) double(solxyz2.y(1,1))...
%     double(solxyz2.z(1,1))];
% fp2q7= [double(solxyz2.x(2,1)) double(solxyz2.y(2,1))...
%     double(solxyz2.z(2,1))];
% fp3q7= [double(solxyz2.x(3,1)) double(solxyz2.y(3,1))...
%     double(solxyz2.z(3,1))];
% 
% return
% %%
% func = @(t,x) [(x(1) + x(2) - q*x(1).^2 - x(1)*x(2) ) / e; -x(2) + f*x(3) - x(1)*x(2) ; (x(1)-x(3)) / p];
% % [t,xa] = ode45(func,[0 100],[0.5 0.5 0.5]);
% [t,xa] = ode45(func,[0 .10],[1400    1.5    1400]);
% 
% plot3(xa(:,1),xa(:,2),xa(:,3));
% hold on
% 
% var1 =[ 0 0 0];
% var2 =[ 1413.71          1.00       1413.71];
% 
% % 
% % var2 =[ 19/2 - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f,...
% %     (3*f)/4 + (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80, 19/2 ...
% %     - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f];
% % var3 =[ (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2, (3*f)/4 ...
% %     - (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80,...
% %     (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2];
% 
% 
% plot3(var1(1),var1(2),var1(3),'r*')
% plot3(var2(1),var2(2),var2(3),'g*')
% 
% % plot3(var1(1),var1(2),var1(3),'r*')
% % plot3(var2(1),var2(2),var2(3),'g*')
% 
% %plot3(var3(1),var3(2),var3(3),'b*')
% 
% grid on
% title('Solution curve')
% %%
% 
% 
% % attempt on quiver
% x= -3:0.5:3;
% y= -3:0.5:3;
% z= -3:0.5:3;
% [X,Y,Z] = meshgrid(x, y, z);
% e= 10^-4; f= 1; p= 100; q= 10^-6; % condition 2
% u=0 ; % without external excitation
% eq11= (x + y - q.*x.^2 - x.*y + u) / e == 0;% 1st order system x
% eq12= -y + f.*z - x.*y + u == 0; % 1st order system y
% eq13= (x-z) / p == 0; % 1st order system z
% % using simplified equation
% h= 100.*x - 5.*x.^2 + (100.*f.*x)/(x + 1) - (100.*f.*x.^2)/(x + 1);
% i= 100.*x - 5.*x.^2 + (100.*f.*x)/(x + 1) - (100.*f.*x.^2)/(x + 1);
% j= 100.*x - 5.*x.^2 + (100.*f.*x)/(x + 1) - (100.*f.*x.^2)/(x + 1);
% [R,S,T] = meshgrid(eq11, eq12, eq13);
% [H,I,J] = meshgrid(h, i, j);
% figure
% quiver3(X,Y,Z,R,S,T)
% view(-35,45)
% figure
% quiver3(X,Y,Z,H,I,J)
% view(-35,45)

%% QUESTION 8 IMPLEMENT CONDITION 3
% close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; f= 1.2; p= 0.5; q= 0.05; % condition 3
u=0 ; % without external excitation
func = @(t,x) [(x(1) + x(2) - q*x(1).^2 - x(1)*x(2) ) / e; -x(2) + f*x(3) - x(1)*x(2) ; (x(1)-x(3)) / p];
% initial condition 1
[t,xa] = ode45(func,[0 100],[20 0.6 6]);
plot3(xa(:,1),xa(:,2),xa(:,3)); hold on
% initial condition 2
[t,xa] = ode45(func,[0 100],[3.8 7.8 6.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
% initial condition 3
[t,xa] = ode45(func,[0 100],[0.5 0.5 0.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
% initial condition 4
[t,xa] = ode45(func,[0 100],[5 4.5 2.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
% initial condition 5
[t,xa] = ode45(func,[0 100],[3.5 1 2.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
hold on
var1 =[ 0 0 0];
var2 =[ 19/2 - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f,...
    (3*f)/4 + (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80, 19/2 ...
    - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f];
var3 =[ (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2, (3*f)/4 ...
    - (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80,...
    (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2];

plot3(var1(1),var1(2),var1(3),'r*')
plot3(var2(1),var2(2),var2(3),'g*')
plot3(var3(1),var3(2),var3(3),'b*')
plot3(20, 0.6, 6,'k*')
plot3(3.8, 7.8, 6.5,'k*')
plot3(0.5, 0.5, 0.5,'k*')
plot3(5, 4.5, 2.5,'k*')
plot3(3.5, 1, 2.5,'k*')
grid on
title('Plot for f=1,2')
%% QUESTION 9 IMPLEMENT DIFFERENT f
% close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; p= 0.5; q= 0.05; % condition 3
f= 0.1; % different value of f
% f= 1.4;
u=0 ; % without external excitation
func = @(t,x) [(x(1) + x(2) - q*x(1).^2 - x(1)*x(2) ) / e; -x(2) + f*x(3) - x(1)*x(2) ; (x(1)-x(3)) / p];
% initial condition 1
[t,xa] = ode45(func,[0 100],[20 0.6 6]);
plot3(xa(:,1),xa(:,2),xa(:,3)); hold on
% initial condition 2
[t,xa] = ode45(func,[0 100],[3.8 7.8 6.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
% initial condition 3
[t,xa] = ode45(func,[0 100],[0.5 0.5 0.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
% initial condition 4
[t,xa] = ode45(func,[0 100],[5 4.5 2.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
% initial condition 5
[t,xa] = ode45(func,[0 100],[3.5 1 2.5]);
plot3(xa(:,1),xa(:,2),xa(:,3));
hold on
var1 =[ 0 0 0];
var2 =[ 19/2 - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f,...
    (3*f)/4 + (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80, 19/2 ...
    - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f];
var3 =[ (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2, (3*f)/4 ...
    - (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80,...
    (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2];

plot3(var1(1),var1(2),var1(3),'r*')
plot3(var2(1),var2(2),var2(3),'g*')
plot3(var3(1),var3(2),var3(3),'b*')
plot3(20, 0.6, 6,'k*')
plot3(3.8, 7.8, 6.5,'k*')
plot3(0.5, 0.5, 0.5,'k*')
plot3(5, 4.5, 2.5,'k*')
plot3(3.5, 1, 2.5,'k*')
grid on
title('Plot for f=1,2')
%% QUESTION 10
close all, clear all, clc;
syms e f p q u x y z xs ys a
e= 10^-2; p= 0.5; q=0.05; f= 1.2;% condition 3
% x= 0; y=0; z=0; % fixed point 1
% x= -9.589; y= 1.34; z= -9.589;% fixed point 2
x= 4.589; y= 0.985; z= 4.589; % fixed point 3
ux= a*(x-xs); uy= a*(y-ys); % with external excitation
eq101= (x + y - q*x^2 - x*y + ux) / e == 0;% 3rd order system x
eq102= -y + f*z - x*y + uy == 0; % 3rd order system y
eq103= (x-z) / p == 0; % 3rd order system z
am= -0.05:0.01:0.05; % values of alpha
solxyz3= solve([(x + y - q*x^2 - x*y + ux) / e,...
    -y + f*z - x*y + uy, (x-z) / p], [xs,ys]); % solve for xs and ys
for i= 1:numel(am)
    a= am(i);
    if a == 0
        xsr(i)= x;
        ysr(i)= y;
        zsr(i)= z;
    else
    xsr(i)= double(subs(solxyz3.xs));
    ysr(i)= double(subs(solxyz3.ys));
    end
end
%% validating
close all, clear all, clc;
syms e f p q u x y z xs ys a
ux= a*(x-xs); uy= a*(y-ys); % with external excitation
eq101= (x + y - q*x^2 - x*y + ux) / e;% 3rd order system x
eq102= -y + f*z - x*y + uy; % 3rd order system y
eq103= (x-z) / p; % 3rd order system z
variable=jacobian([eq101, eq102, eq103],[x, y, z]);
% x= 0; y=0; z=0; % fixed point 1
% x= -9.589; y= 1.34; z= -9.589;% fixed point 2
x= 4.589; y= 0.985; z= 4.589; % fixed point 3
e= 10^-2; p= 0.5; q=0.05; f= 1.2;% condition 3
am= -0.05:0.001:0.05; % values of alpha
count=0;
for i= 1:numel(am)
    a= am(i);
    x11= double(subs(variable(1,1)));
    x12= double(subs(variable(1,2)));
    x13= double(subs(variable(1,3)));
    x21= double(subs(variable(2,1)));
    x22= double(subs(variable(2,2)));
    x23= double(subs(variable(2,3)));
    x31= double(subs(variable(3,1)));
    x32= double(subs(variable(3,2)));
    x33= double(subs(variable(3,3)));
    mat{i}= [x11 x12 x13; x21 x22 x23; x31 x32 x33];
    ev{i}= eig(mat{i});
    tev(i,1)= ev{i}(1);
    tev(i,2)= ev{i}(2);
    tev(i,3)= ev{i}(3);
    if ev{i}(2)<0
        count=count+1;
        am(count)
    end
end