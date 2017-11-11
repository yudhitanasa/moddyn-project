%=========================================================================
% I.Y. Tanasa (1034117)       
% Aleman Zapata, R.A. (59383)
% Modeling Dynamics Project
%=========================================================================
%% QUESTION 8 IMPLEMENT CONDITION 3
% close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; f= 1.2; p= 0.5; q= 0.05; % condition 3
u=0 ; % without external excitation
func = @(t,x) [(x(1) + x(2) - q*x(1).^2 - x(1)*x(2) ) / e; -x(2) + f*x(3) - x(1)*x(2) ; (x(1)-x(3)) / p];
[t,xa] = ode45(func,[0 100],[0.5 0.5 0.5]);
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

grid on
title('Plot for f=1,2')
