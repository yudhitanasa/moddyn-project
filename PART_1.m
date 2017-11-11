%=========================================================================
% I.Y. Tanasa (1034117)       
% Aleman Zapata, R.A. (59383)
% Modeling Dynamics Project for question 1 to 3
%=========================================================================
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
fp1q1= [solxyz1.x(1,1) solxyz1.y(1,1) solxyz1.z(1,1)]; % (0,0,0)
fp2q1= [solxyz1.x(2,1) solxyz1.y(2,1) solxyz1.z(2,1)]; % depend on f
fp3q1= [solxyz1.x(3,1) solxyz1.y(3,1) solxyz1.z(3,1)]; % depend on f
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
    fp1q1f(i,1)= double(subs(solxyz1.x(1,1))); % table fp 1
    fp1q1f(i,2)= double(subs(solxyz1.y(1,1))); 
    fp1q1f(i,3)= double(subs(solxyz1.z(1,1)));
    fp2q1f(i,1)= double(subs(solxyz1.x(2,1))); % table fp 2 
    fp2q1f(i,2)= double(subs(solxyz1.y(2,1))); 
    fp2q1f(i,3)= double(subs(solxyz1.z(2,1)));
    fp3q1f(i,1)= double(subs(solxyz1.x(3,1))); % table fp 3 
    fp3q1f(i,2)= double(subs(solxyz1.y(3,1)));
    fp3q1f(i,3)= double(subs(solxyz1.z(3,1)));
end
ftab1q1= horzcat((fm'),fp1q1f); % final table with f,y,z
ftab2q1= horzcat((fm'),fp2q1f);
ftab3q1= horzcat((fm'),fp3q1f);
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
ftablam1q2= horzcat((fm'),lamfp1q2); % final table with f,x,y,z
ftablam2q2= horzcat((fm'),lamfp2q2);
ftablam3q2= horzcat((fm'),lamfp3q2);
format bank
retab1= horzcat((fm'),fp1q1f,lamfp1q2);
retab2= horzcat((fm'),fp2q1f,lamfp2q2);
retab3= horzcat((fm'),fp3q1f,lamfp3q2);
%% QUESTION 3 FIND BIFURCATION POINTS
% analyticaly using table we found at question 2
