close all, clear all, clc;
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


%%
func = @(t,x) [(x(1) + x(2) - q*x(1).^2 - x(1)*x(2) ) / e; -x(2) + f*x(3) - x(1)*x(2) ; (x(1)-x(3)) / p];
% [t,xa] = ode45(func,[0 100],[0.5 0.5 0.5]);
[t,xa1] = ode45(func,[0 0.5],[4    1.5    1400]);
[t,xa2] = ode45(func,[0 0.5],[5    1    2000]);
[t,xa3] = ode45(func,[0 0.5],[0.1    600    100]);
[t,xa4] = ode45(func,[0 0.5],[0.01    100    600]);

%%
figure()
plot(xa1(:,3))
hold on
plot(xa2(:,3))
plot(xa3(:,3))
plot(xa4(:,3))

%%

figure()
plot(xa1(:,2))
hold on
plot(xa2(:,2))
plot(xa3(:,2))
plot(xa4(:,2))

%%

figure()
plot(xa1(:,1))
hold on
plot(xa2(:,1))
plot(xa3(:,1))
plot(xa4(:,1))
%% FInding Y 
[t,xa4] = ode15s(func,[0 0.5],[0.1    100    100]);
[t,xa5] = ode15s(func,[0 0.5],[0.5    200    200]);
[t,xa6] = ode15s(func,[0 0.5],[1    300    300]);
[t,xa7] = ode15s(func,[0 0.5],[0.2    400    400]);
[t,xa8] = ode15s(func,[0 0.5],[1.5    500    500]);
[t,xa9] = ode15s(func,[0 0.5],[2    600    600]);
%%


figure()
plot(xa4(:,2))
hold on
plot(xa5(:,2))
plot(xa6(:,2))
plot(xa7(:,2))
plot(xa8(:,2))
plot(xa9(:,2))

title('Phase Plot for Y')
xlabel('Time (sec)')
ylabel('Magnitude')
%%

figure()
plot(xa4(:,1))
hold on
plot(xa5(:,1))
plot(xa6(:,1))
plot(xa7(:,1))
plot(xa8(:,1))
plot(xa9(:,1))

title('Phase Plot for X')
xlabel('Time (sec)')
ylabel('Magnitude')

%%
%%

figure()
plot(xa4(:,3))
hold on
plot(xa5(:,3))
plot(xa6(:,3))
plot(xa7(:,3))
plot(xa8(:,3))
plot(xa9(:,3))

title('Phase Plot for Z')
xlabel('Time (sec)')
ylabel('Magnitude')
%%
plot3(xa4(:,1),xa4(:,2),xa4(:,3));
hold on
plot3(xa5(:,1),xa5(:,2),xa5(:,3));
plot3(xa6(:,1),xa6(:,2),xa6(:,3));
plot3(xa7(:,1),xa7(:,2),xa7(:,3));
plot3(xa8(:,1),xa8(:,2),xa8(:,3));
plot3(xa9(:,1),xa9(:,2),xa9(:,3));
%%
plot3(xa1(:,1),xa1(:,2),xa1(:,3));
hold on
plot3(xa2(:,1),xa2(:,2),xa2(:,3));
plot3(xa3(:,1),xa3(:,2),xa3(:,3));



% [t,xa1] = ode45(func,[0 .1],[400    1.5    1400]);
% [t,xa2] = ode45(func,[0 .1],[500    1    2000]);
% [t,xa3] = ode45(func,[0 .1],[0.1    300    100]);

plot3(400,1.5 ,1400,'r*')
plot3(500,1,2000,'g*')
plot3(0.1,300,100,'b*')
% var1 =[ 0 0 0];
% var2 =[ 1413.71          1.00       1413.71];
% 
% plot3(var1(1),var1(2),var1(3),'r*')
% plot3(var2(1),var2(2),var2(3),'g*')

%%

% 
% var2 =[ 19/2 - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f,...
%     (3*f)/4 + (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80, 19/2 ...
%     - (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f];
% var3 =[ (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2, (3*f)/4 ...
%     - (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80,...
%     (400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2];



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
