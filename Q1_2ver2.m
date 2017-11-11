close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
% eq11= (x + y - q*x^2 - x*y + u) / e == 0;% 1st order system x
eq11= x + y - q*x^2 - x*y + u == 0;% 1st order system x
eq12= -y + f*z - x*y + u == 0; % 1st order system y
eq13= (x-z) / p == 0; % 1st order system z
solz1= solve(eq13, z);
eq12ex= -y + f*solz1 - x*y + u == 0; % eliminate z
soly2= solve(eq12ex, y);
eq11ex= (x + soly2 - q*x^2 - x*soly2 + u) / e == 0; %eliminate y
solx2= solve(eq11ex, x);
for x=solx2(1)
    eq13x1= (x-z) / p == 0;
    solz1x1= solve(eq13x1, z);
    eq12x1= -y + f*solz1x1 - x*y + u == 0;
    soly1x1= solve(eq12x1, y);
    fp1q1=[solx2(1) soly1x1 solz1x1];
    % fixed point set 1 (x,y,z) independent of f with x= 0
end
for x=solx2(2)
    eq13x2= (x-z) / p == 0;
    solz1x2= solve(eq13x2, z);
    eq12x2= -y + f*solz1x2 - x*y + u == 0;
    soly1x2= solve(eq12x2, y);
    fp2q1=[solx2(2) soly1x2 solz1x2];
    % fixed point set 2 (x,y,z) depend on f with x = 1 of the roots
end
for x=solx2(3)
    eq13x3= (x-z) / p == 0;
    solz1x3= solve(eq13x3, z);
    eq12x3= -y + f*solz1x2 - x*y + u == 0;
    soly1x3= solve(eq12x3, y);
    fp3q1=[solx2(3) soly1x3 solz1x3];
    % fixed point set 3 (x,y,z) depend on f with x = the other roots
end
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
    fp1q1f(i,1)= double(subs(solx2(1))); % table for fp 1 from question 1
    fp1q1f(i,2)= double(subs(soly1x1));
    fp1q1f(i,3)= double(subs(solz1x1));
    fp2q1f(i,1)= double(subs(solx2(2))); % table for fp 2 from question 1
    fp2q1f(i,2)= double(subs(soly1x2)); 
    fp2q1f(i,3)= double(subs(solz1x2));
    fp3q1f(i,1)= double(subs(solx2(3))); % table for fp 3 from question 1
    fp3q1f(i,2)= double(subs(soly1x3));
    fp3q1f(i,3)= double(subs(solz1x3));
end
ftab1q1= horzcat((fm'),fp1q1f); % final table with f,x,y,z
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
        x21q2(i,1), x22q2(i,1), x23q2(i,1);...
        x31q2(i,1), x32q2(i,1), x33q2(i,1)]; % each cell contain jacobian matrix
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
        x21q2(i,1), x22q2(i,1), x23q2(i,1);...
        x31q2(i,1), x32q2(i,1), x33q2(i,1)]; % each cell contain jacobian matrix
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
        x21q2(i,1), x22q2(i,1), x23q2(i,1);...
        x31q2(i,1), x32q2(i,1), x33q2(i,1)]; % each cell contain jacobian matrix
end
for j=1:numel(fm)
    eigfp1q2{j,1}= eig(mfp1q2{j,1}); % each cell contain lambda 1,lambda 2 
    eigfp2q2{j,1}= eig(mfp2q2{j,1}); % and lamda 3 
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
















