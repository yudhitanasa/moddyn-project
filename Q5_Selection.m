close all, clear all, clc;
syms e f p q u x y z  % redefine
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation

% USING X = NEGATIVE VALUE (RESULT FROM solx(1)------equation 9)
eq22m= -y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))*y == 0; % 2nd order system y dot, with z=x
eq23m= ((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p == 0; % 2nd order system z dot
solyzm= solve([eq22m, eq23m ], [y,z]) % set of fixed point (y,z)

% USING X = POSITIVE VALUE (RESULT FROM solx(2)------equation 10)
eq22= -y+f*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))*y == 0; % 2nd order system y dot, with z=x
eq23= ((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))-z)/p == 0; % 2nd order system z dot
solyz= solve([eq22, eq23 ], [y,z]); % set of fixed point (y,z)

e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
clc;
A=jacobian([eq22 eq23],[y,z]);
B=jacobian([eq22m eq23m],[y,z]);
pretty(A(2,1))
pretty(B(2,1))

return
f=0.1; % because fp are same we only test the one from negative x
fp1q5y= double(subs(solyzm.y(1,1)))
fp1q5z= double(subs(solyzm.z(1,1)))    
fp2q5y= double(subs(solyzm.y(2,1)))
fp2q5z= double(subs(solyzm.z(2,1)))
fp3q5y= double(subs(solyzm.y(3,1)))
fp3q5z= double(subs(solyzm.z(3,1)))
% FP1
% y= fp1q5y;
% z= fp1q5z;
% FP2
% y= fp2q5y;
% z= fp2q5z;
% % FP3
y= fp3q5y;
z= fp3q5z;
double(subs(-y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))*y))
double(subs(((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p))

















