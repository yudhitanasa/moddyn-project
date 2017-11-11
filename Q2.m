%=========================================================================
% I.Y. Tanasa (1034117)       
% Aleman Zapata, R.A. (59383)
% Modeling Dynamics Project
%=========================================================================
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
