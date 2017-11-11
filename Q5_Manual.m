close all, clear all, clc;
syms e f p q u x y z  % redefine
% e= 10^-2; p= 0.5; q=0.05; % condition 1
% u=0 ; % without external excitation
eq22= -y+f*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))*y == 0; % 2nd order system y dot, with z=x
eq23= ((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))-z)/p == 0; % 2nd order system z dot
solyz= solve([eq22, eq23 ], [y,z]); % two set of fixed point (y,z)
% e= 10^-2; p= 0.5; q=0.05; % condition 1
% u=0 ; % without external excitation
eq22m= -y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))*y == 0; % 2nd order system y dot, with z=x
eq23m= ((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p == 0; % 2nd order system z dot
solyzm= solve([eq22m, eq23m ], [y,z]) % two set of fixed point (y,z)
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
fm= (0.1:0.1:2);  % condition 1
% jacobian([eq22 eq23],[y,z])
jacobian([eq22m eq23m],[y,z])
% return
for i=1:numel(fm)
    y= fp1q5f(i,1); % for 1st fp from question 5
    x11(i,1)= (y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q)...
        + (y*((4*q + 2*y - 2)/(2*(4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2))...
        + 1))/(2*q) - 1;
    x12(i,1)= fm(i);
    x21(i,1)=  -((4*q + 2*y - 2)/(2*(4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2))...
        + 1)/(2*p*q) ;
    x22(i,1)= -1/p;
    mfp1{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
%     ====================================================================    
    y= fp2q5f(i,1); % for 2nd fp from question 5
    x11(i,2)= (y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q)...
        + (y*((4*q + 2*y - 2)/(2*(4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2))...
        + 1))/(2*q) - 1;
    x12(i,2)= fm(i);
    x21(i,2)= -((4*q + 2*y - 2)/(2*(4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2))...
        + 1)/(2*p*q) ;
    x22(i,2)= -1/p;
    mfp2{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
end