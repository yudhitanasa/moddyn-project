%% QUESTION 1 FIND ROOT OF 3 VARIABLES EXPRESSION x, y AND z
close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
eq11= (x + y - q*x^2 - x*y + u) / e == 0;% 1st order system x
% eq11= x + y - q*x^2 - x*y + u == 0;% 1st order system x
eq12= -y + f*z - x*y + u == 0; % 1st order system y
eq13= (x-z) / p == 0; % 1st order system z
solxyz1= solve([(x + y - q*x^2 - x*y + u) / e,...
    -y + f*z - x*y + u, (x-z) / p], [x,y,z]); % 3 set of fp (x,y,z)
% solxyz1= solve([eq11, eq12, eq13 ], [x,y,z]); % 3 set of fp (x,y,z)
fp1q1= [solxyz1.x(1,1) solxyz1.y(1,1) solxyz1.z(1,1)]; % depend on f
fp2q1= [solxyz1.x(2,1) solxyz1.y(2,1) solxyz1.z(2,1)]; % depend on f
fp3q1= [solxyz1.x(3,1) solxyz1.y(3,1) solxyz1.z(3,1)]; % (0,0)
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
    fp1q1f(i,1)= double(subs(solxyz1.x(1,1))); % table fp 1 from question 4
    fp1q1f(i,2)= double(subs(solxyz1.y(1,1))); 
    fp1q1f(i,3)= double(subs(solxyz1.z(1,1)));
    fp2q1f(i,1)= double(subs(solxyz1.x(2,1))); % table fp 2 from question 4
    fp2q1f(i,2)= double(subs(solxyz1.y(2,1))); 
    fp2q1f(i,3)= double(subs(solxyz1.z(2,1)));
    fp3q1f(i,1)= double(subs(solxyz1.x(3,1))); % table fp 3 from question 4
    fp3q1f(i,2)= double(subs(solxyz1.y(3,1)));
    fp3q1f(i,3)= double(subs(solxyz1.z(3,1)));
end
ftab1q1= horzcat((fm'),fp1q1f); % final table with f,y,z
ftab2q1= horzcat((fm'),fp2q1f);
ftab3q1= horzcat((fm'),fp3q1f);
%{
%% QUESTION 2 ASSES STABILITY OF THE FIXED POINTS
% close all, clear all, clc;
syms e f p q u x y z % redefine
% we use eq11, eq12 and eq13 from before
jacobian([(x + y - q*x^2 - x*y + u) / e, -y + f*z - x*y + u,...
    (x-z) / p],[x,y,z]); % check jacobian here
% jacobian([eq11, eq12, eq13],[x,y,z]) % check jacobian here
e= 10^-2; p= 0.5; q=0.05; % condition 1
fm= (0.1:0.1:2);  % condition 1
% jacobian matrix plug in 3 fp and every f
for i=1:numel(fm)
    x= fp1q1f(i,1); % x value for 1st fp from question 1
    y= fp1q1f(i,2); % y value for 1st fp from question 1
    x11q2(i,1)= -(y + 2*q*x - 1)/e;
    x12q2(i,1)= -(x - 1)/e;
    x13q2(i,1)= 0;
    x21q2(i,1)= -y;
    x22q2(i,1)= - x - 1;
    x23q2(i,1)= fm(i);
    x31q2(i,1)= 1/p;
    x32q2(i,1)= 0;
    x33q2(i,1)= -1/p;
    mfp1q2{i,1}= [x11q2(i,1), x12q2(i,1), x13q2(i,1);...
        x21q2(i,1), x22q2(i,1), x23q2(i,1);... % each cell contain
        x31q2(i,1), x32q2(i,1), x33q2(i,1)];   % jacobian matrix
%     ===================================================================
    x= fp2q1f(i,1); % x value for 2nd fp from question 1
    y= fp2q1f(i,2); % y value for 2nd fp from question 1
    x11q2(i,1)= -(y + 2*q*x - 1)/e;
    x12q2(i,1)= -(x - 1)/e;
    x13q2(i,1)= 0;
    x21q2(i,1)= -y;
    x22q2(i,1)= - x - 1;
    x23q2(i,1)= fm(i);
    x31q2(i,1)= 1/p;
    x32q2(i,1)= 0;
    x33q2(i,1)= -1/p;
    mfp2q2{i,1}= [x11q2(i,1), x12q2(i,1), x13q2(i,1);...
        x21q2(i,1), x22q2(i,1), x23q2(i,1);... % each cell contain
        x31q2(i,1), x32q2(i,1), x33q2(i,1)];   % jacobian matrix
%     ===================================================================
    x= fp3q1f(i,1); % x value for 3rd fp from question 1
    y= fp3q1f(i,2); % y value for 3rd fp from question 1
    x11q2(i,1)= -(y + 2*q*x - 1)/e;
    x12q2(i,1)= -(x - 1)/e;
    x13q2(i,1)= 0;
    x21q2(i,1)= -y;
    x22q2(i,1)= - x - 1;
    x23q2(i,1)= fm(i);
    x31q2(i,1)= 1/p;
    x32q2(i,1)= 0;
    x33q2(i,1)= -1/p;
    mfp3q2{i,1}= [x11q2(i,1), x12q2(i,1), x13q2(i,1);...
        x21q2(i,1), x22q2(i,1), x23q2(i,1);... % each cell contain
        x31q2(i,1), x32q2(i,1), x33q2(i,1)];   % jacobian matrix
end
for j=1:numel(fm)
    eigfp1q2{j,1}= eig(mfp1q2{j,1}); % each cell contain lambda 1, 
    eigfp2q2{j,1}= eig(mfp2q2{j,1}); % lambda 2 and lamda 3 
    eigfp3q2{j,1}= eig(mfp3q2{j,1});
%     ================================================================
    lamfp1q2(j,1)= eigfp1q2{j,1}(1);  % for simplification
    lamfp1q2(j,2)= eigfp1q2{j,1}(2);
    lamfp1q2(j,3)= eigfp1q2{j,1}(3);
    lamfp2q2(j,1)= eigfp2q2{j,1}(1);
    lamfp2q2(j,2)= eigfp2q2{j,1}(2);
    lamfp2q2(j,3)= eigfp2q2{j,1}(3);
    lamfp3q2(j,1)= eigfp3q2{j,1}(1);
    lamfp3q2(j,2)= eigfp3q2{j,1}(2);
    lamfp3q2(j,3)= eigfp3q2{j,1}(3);
end
ftab1q2= horzcat((fm'),lamfp1q2); % final table with f,x,y,z
ftab2q2= horzcat((fm'),lamfp2q2);
ftab3q2= horzcat((fm'),lamfp3q2);
%% QUESTION 3 FIND BIFURCATION POINTS
% analyticaly using table we found at question 2
%}
%% QUESTION 4 DERIVE 2ND ORDER DYNAMICAL SYSTEM
% close all, clear all, clc;
syms e f p q u x y z  % redefine
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
% we want to eliminate eq11, we use define x as
solx= solve(eq11,x)
return
sol=2; % i still dont know why Daan use (i=1) we have 2 of solx
eq23= (solx(sol)-z)/p == 0; % 2nd order system z dot
solz= solve(eq23, z); % we define z to be used in eq22 (p is not zero)
eq22= -y+f*solz-solx(sol)*y == 0; % 2nd order system y dot, with z=x
% soly= solve(eq22, y) % 2nd order system y dot

jacobian([eq22 eq23],[y,z]) % check jacobian here



%% QUESTION 5 FIND ROOT OF 2 VARIABLES EXPRESSION y AND z
solyz= solve([y*(10*y + 2*5^(1/2)*(5*y^2 - 9*y + 5)^(1/2) - 10)...
    - f*(10*y + 2*5^(1/2)*(5*y^2 - 9*y + 5)^(1/2) - 10) - y,...
    -(10*y + z + 2*5^(1/2)*(5*y^2 - 9*y + 5)^(1/2) - 10)/p],...
    [y,z]); % two set of fixed point (y,z)
% solyz= solve([eq22, eq23 ], [y,z]); % two set of fixed point (y,z)
fp1q5= [solyz.y(1,1) solyz.z(1,1)]; % depend on f
fp2q5= [solyz.y(2,1) solyz.z(2,1)]; % depend on f
fp3q5= [solyz.y(3,1) solyz.z(3,1)]; % (0,0)
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
    fp1q5f(i,1)= double(subs(solyz.y(1,1))); % table fp 1 from question 4
    fp1q5f(i,2)= double(subs(solyz.z(1,1)));    
    fp2q5f(i,1)= double(subs(solyz.y(2,1))); % table fp 2 from question 4
    fp2q5f(i,2)= double(subs(solyz.z(2,1)));    
    fp3q5f(i,1)= double(subs(solyz.y(3,1))); % table fp 3 from question 4
    fp3q5f(i,2)= double(subs(solyz.z(3,1)));
end
tab1q5= horzcat((fm'),fp1q5f); % table with f,y,z
tab2q5= horzcat((fm'),fp2q5f);
tab3q5= horzcat((fm'),fp3q5f);
%{
%% QUESTION 5 FIND JACOBIAN TO ASSESS STABILITY
% close all, clear all, clc;
syms e f p q u x y z % redefine
% we use eq22 and eq23 from before
% jacobian([y*(10*y + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 10)...
%     - f*(10*y + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 10) - y,...     % check
%     -(10*y + z + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 10)/p],[y,z]); % here
jacobian([eq22 eq23],[y,z]) % check jacobian here
p= 0.5; fm= (0.1:0.1:2);  % condition 1
% jacobian matrix plug in 3 values of y i.e. from fp 1 to fp3 
% and 20 values of f 
for i=1:numel(fm)
    y= fp1q5f(i,1); % for 1st fp from question 5
    x11(i,1)= 10*y - fm(i)*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)...
        + y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)...
        + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
    x12(i,1)= 0;
    x21(i,1)= -((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)/p;
    x22(i,1)= -1/p;
    mfp1{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
%     ====================================================================    
    y= fp2q5f(i,1); % for 2nd fp from question 5
    x11(i,2)= 10*y - fm(i)*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)...
        + y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)...
        + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
    x12(i,2)= 0;
    x21(i,2)= -((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)/p;
    x22(i,2)= -1/p;
    mfp2{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
%     ====================================================================    
    y= fp3q5f(i,1); % for 3rd fp from question 5
    x11(i,3)= 10*y - fm(i)*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)...
        + y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)...
        + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
    x12(i,3)= 0;
    x21(i,3)= -((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)/p;
    x22(i,3)= -1/p;
    mfp3{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
end
% find the eigen value
for j=1:numel(fm)
eigfp1{j,1}= eig(mfp1{j,1}); % each cell contain lambda 1
eigfp2{j,1}= eig(mfp2{j,1}); % and lambda 2 for 1 value f
eigfp3{j,1}= eig(mfp3{j,1});
lam1fp1(j,1)= eigfp1{j,1}(1);
lam2fp1(j,2)= eigfp1{j,1}(2);
lam1fp2(j,1)= eigfp2{j,1}(1);
lam2fp2(j,2)= eigfp2{j,1}(2);
lam1fp3(j,1)= eigfp3{j,1}(1);
lam2fp3(j,2)= eigfp3{j,1}(2);
end
ftab1q5= horzcat(tab1q5,lam1fp1,lam2fp1); % final table with f,y,z,
ftab2q5= horzcat(tab2q5,lam1fp2,lam2fp2); % lambda 1, lambda 2
ftab3q5= horzcat(tab3q5,lam1fp3,lam2fp3);

%% QUESTION 6 FIND BIFURCATION POINTS
% analyticaly using table we found at question 5
return
%% QUESTION 7 VERIVY THE POSITIVE INVARIANT FOR UNDRIVEN SYSTEM
% close all, clear all, clc;
syms e f p q u x y z
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
%}