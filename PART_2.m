%=========================================================================
% I.Y. Tanasa (1034117)       
% Aleman Zapata, R.A. (59383)
% Modeling Dynamics Project for question 4 to 6
%=========================================================================
%% QUESTION 4 DERIVE 2ND ORDER DYNAMICAL SYSTEM
close all, clear all, clc;
syms e f p q u x y z  % redefine
% e= 10^-2; p= 0.5; q=0.05; % condition 1
% u=0 ; % without external excitation
eq11= (x + y - q*x^2 - x*y + u) / e == 0;% 1st order system x
% we want to eliminate eq11, we use define x as
solx= solve(eq11,x);
sol=2; 
eq23= (solx(sol)-z)/p == 0; % 2nd order system z dot
solz= solve(eq23, z); % we define z to be used in eq22 (p is not zero)
eq22= -y+f*solz-solx(sol)*y == 0; % 2nd order system y dot, with z=x
% soly= solve(eq22, y) % 2nd order system y dot
%% QUESTION 5 FIND ROOT OF 2 VARIABLES EXPRESSION y AND z
% close all, clear all, clc;
syms e f p q u x y z  % redefine
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
% USING X = NEGATIVE VALUE (RESULT FROM solx(1)------equation 9)
% 2nd order system y dot, with z=x
eq22n= -y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2)-1)/(2*q))*y == 0;
% 2nd order system z dot
eq23n= ((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p == 0;
solyzn= solve([eq22n, eq23n ], [y,z]); % set of fixed point (y,z)
fp1nq5= [solyzn.y(1,1) solyzn.z(1,1)]; % (0,0)
fp2nq5= [solyzn.y(2,1) solyzn.z(2,1)]; % depend on f
% fp3nq5= [solyz.y(3,1) solyzn.z(3,1)]; % depend on f BUT INVALID
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
    fp1nq5f(i,1)= double(subs(solyzn.y(1,1))); % table fp 1 from question 4
    fp1nq5f(i,2)= double(subs(solyzn.z(1,1)));    
    fp2nq5f(i,1)= double(subs(solyzn.y(2,1))); % table fp 2 from question 4
    fp2nq5f(i,2)= double(subs(solyzn.z(2,1)));    
%     fp3nq5f(i,1)= double(subs(solyzn.y(3,1))); % INVALID NOT VERIFIED
%     fp3nq5f(i,2)= double(subs(solyzn.z(3,1))); % INVALID NOT VERIFIED
end
tab1nq5= horzcat((fm'),fp1nq5f); % table with f,y,z
tab2nq5= horzcat((fm'),fp2nq5f);
% tab3nq5= horzcat((fm'),fp3nq5f); % INVALID NOT VERIFIED
% ========================================================================
syms e f p q u x y z  % redefine
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
fm= (0.1:0.1:2);  % condition 1
% USING X = POSITIVE VALUE (RESULT FROM solx(2)------equation 10)
% 2nd order system y dot, with z=x
eq22p= -y+f*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2)-y + 1)/(2*q))*y == 0;
% 2nd order system z dot
eq23p= ((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2)-y + 1)/(2*q))-z)/p == 0;
solyzp= solve([eq22p, eq23p], [y,z]); % set of fixed point (y,z)
% fp1pq5= [solyzp.y(1,1) solyzp.z(1,1)]; % n/a
fp2pq5= [solyzp.y(1,1) solyzp.z(1,1)]; % depend on f
fp3pq5= [solyzp.y(2,1) solyzp.z(2,1)]; % depend on f BUT INVALID
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
%     fp1pq5f(i,1)= double(subs(solyzp.y(1,1))); % n/a
%     fp1pq5f(i,2)= double(subs(solyzp.z(1,1))); % n/a
    fp2pq5f(i,1)= double(subs(solyzp.y(1,1))); % table fp 2 from question 4
    fp2pq5f(i,2)= double(subs(solyzp.z(1,1)));    
    fp3pq5f(i,1)= double(subs(solyzp.y(2,1))); % INVALID NOT VERIFIED
    fp3pq5f(i,2)= double(subs(solyzp.z(2,1))); % INVALID NOT VERIFIED
end
% tab1pq5= horzcat((fm'),fp1pq5f); % n/a
tab2pq5= horzcat((fm'),fp2pq5f);
tab3pq5= horzcat((fm'),fp3pq5f); % INVALID NOT VERIFIED
%% QUESTION 5 FIND JACOBIAN TO ASSESS STABILITY
% close all, clear all, 
clc;
% FOR NEGATIVE SIGN OF SQUARE ROOOT =====================================
syms e f p q u x y z % redefine
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
fm= (0.1:0.1:2);  % condition 1
% USING X = NEGATIVE VALUE (RESULT FROM solx(1)------equation 9)
% 2nd order system y dot, with z=x
eq22n= -y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))*y == 0;
% 2nd order system z dot
eq23n= ((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p == 0;
jacobian([eq22n eq23n],[y,z]); % check jacobian here
% e= 10^-2; p= 0.5; q=0.05; % condition 1
% u=0 ; % without external excitation
% fm= (0.1:0.1:2);  % condition 1
% jacobian matrix plug in 3 values of y i.e. from fp 1 to fp3 
% and 20 values of f 
for i=1:numel(fm)
    y= fp1nq5f(i,1); % for 1st fp from question 5
    x11(i,1)= 10*y + y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10) + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
    x12(i,1)= fm(i);
    x21(i,1)= - (10*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 20;
    x22(i,1)= -1/p;
    mfp1n{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
%     ====================================================================    
    y= fp2nq5f(i,1); % for 2nd fp from question 5
    x11(i,2)= 10*y + y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10) + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
    x12(i,2)= fm(i);
    x21(i,2)= - (10*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 20;
    x22(i,2)= -1/p;
    mfp2n{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
%     ====================================================================    
%     y= fp3nq5f(i,1); % for 3rd fp from question 5
%     x11(i,3)= 10*y + y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) + 10)...
%         + 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
%     x12(i,3)= fm(i);
%     x21(i,3)= - (10*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 20;
%     x22(i,3)= -1/p;
%     mfp3n{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
end
% find the eigen value
for j=1:numel(fm)
eigfp1nq5{j,1}= eig(mfp1n{j,1}); % each cell contain lambda 1
eigfp2nq5{j,1}= eig(mfp2n{j,1}); % and lambda 2 for 1 value f
% eigfp3nq5{j,1}= eig(mfp3n{j,1});
lam1fp1nq5(j,1)= eigfp1nq5{j,1}(1);
lam2fp1nq5(j,1)= eigfp1nq5{j,1}(2);
lam1fp2nq5(j,1)= eigfp2nq5{j,1}(1);
lam2fp2nq5(j,1)= eigfp2nq5{j,1}(2);
% lam1fp3nq5(j,1)= eigfp3nq5{j,1}(1);
% lam2fp3nq5(j,1)= eigfp3nq5{j,1}(2);
end
ftablam1nq5= horzcat(tab1nq5,lam1fp1nq5,lam2fp1nq5); % final table with f,y,z,
ftablam2nq5= horzcat(tab2nq5,lam1fp2nq5,lam2fp2nq5); % lambda 1, lambda 2
% ftablam3nq5= horzcat(tab3nq5,lam1fp3nq5,lam2fp3nq5);


% FOR POSITIVE SIGN OF SQUARE ROOOT =====================================
syms e f p q u x y z % redefine
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
fm= (0.1:0.1:2);  % condition 1
% 2nd order system y dot, with z=x
eq22p= -y+f*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2)-y + 1)/(2*q))*y == 0;
% 2nd order system z dot
eq23p= ((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2)-y + 1)/(2*q))-z)/p == 0;
jacobian([eq22p eq23p],[y,z]); % check jacobian here
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
fm= (0.1:0.1:2);  % condition 1
% jacobian matrix plug in 3 values of y i.e. from fp 1 to fp3 
% and 20 values of f 
% return
for i=1:numel(fm)
%     y= fp2pq5f(i,1); % for 1st fp from question 5
%     x11(i,1)= 10*y - y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 10)...
%         - 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
%     x12(i,1)= fm(i);
%     x21(i,1)= (10*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 20;
%     x22(i,1)= -2;
%     mfp1p{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
%     % each cell contain jacobian matrix
%     ====================================================================    
    y= fp2pq5f(i,1); % for 2nd fp from question 5
    x11(i,2)= 10*y - y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 10) - 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
    x12(i,2)= fm(i);
    x21(i,2)= (10*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 20;
    x22(i,2)= -1/p;
    mfp2p{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
%     ====================================================================    
    y= fp3pq5f(i,1); % for 3rd fp from question 5
    x11(i,3)= 10*y - y*((5*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 10) - 10*(y^2 - (9*y)/5 + 1)^(1/2) - 11;
    x12(i,3)= fm(i);
    x21(i,3)= (10*(2*y - 9/5))/(y^2 - (9*y)/5 + 1)^(1/2) - 20;
    x22(i,3)= -1/p;
    mfp3p{i,1}= [x11(i,1), x12(i,1); x21(i,1), x22(i,1)];
    % each cell contain jacobian matrix
end
% find the eigen value
for j=1:numel(fm)
% eigfp1pq5{j,1}= eig(mfp1p{j,1}); % each cell contain lambda 1
eigfp2pq5{j,1}= eig(mfp2p{j,1}); % and lambda 2 for 1 value f
eigfp3pq5{j,1}= eig(mfp3p{j,1});
% lam1fp1pq5(j,1)= eigfp1pq5{j,1}(1);
% lam2fp1pq5(j,1)= eigfp1pq5{j,1}(2);
lam1fp2pq5(j,1)= eigfp2pq5{j,1}(1);
lam2fp2pq5(j,1)= eigfp2pq5{j,1}(2);
lam1fp3pq5(j,1)= eigfp3pq5{j,1}(1);
lam2fp3pq5(j,1)= eigfp3pq5{j,1}(2);
end
% ftablam1pq5= horzcat(tab1pq5,lam1fp1pq5,lam2fp1pq5); % final table with f,y,z,
ftablam2pq5= horzcat(tab2pq5,lam1fp2pq5,lam2fp2pq5); % lambda 1, lambda 2
ftablam3pq5= horzcat(tab3pq5,lam1fp3pq5,lam2fp3pq5);


%% QUESTION 6 FIND BIFURCATION POINTS
% analyticaly using table we found at question 5
