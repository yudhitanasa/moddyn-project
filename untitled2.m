%% QUESTION 5 FIND ROOT OF 2 VARIABLES EXPRESSION y AND z
close all, clear all, clc;
syms e f p q u x y z  % redefine
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation

% % USING X = NEGATIVE VALUE (RESULT FROM solx(1)------equation 9)
% % 2nd order system y dot, with z=x
% eq22= -y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2)-1)/(2*q))*y == 0;
% % 2nd order system z dot
% eq23= ((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p == 0;

% USING X = POSITIVE VALUE (RESULT FROM solx(2)------equation 10)
% 2nd order system y dot, with z=x
eq22= -y+f*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2)-y + 1)/(2*q))*y == 0;
% 2nd order system z dot
eq23= ((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))-z)/p == 0; 

solyz= solve([eq22, eq23 ], [y,z]); % set of fixed point (y,z)
% fp1q5= [solyz.y(1,1) solyz.z(1,1)]; % (0,0)
% fp2q5= [solyz.y(2,1) solyz.z(2,1)]; % depend on f
fp3q5= [solyz.y(2,1) solyz.z(2,1)]; % depend on f BUT INVALID
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
%     fp1q5f(i,1)= double(subs(solyz.y(1,1))); % table fp 1 from question 4
%     fp1q5f(i,2)= double(subs(solyz.z(1,1)));    
%     fp2q5f(i,1)= double(subs(solyz.y(2,1))); % table fp 2 from question 4
%     fp2q5f(i,2)= double(subs(solyz.z(2,1)));    
    fp3q5f(i,1)= double(subs(solyz.y(2,1))); % INVALID NOT VERIFIED
    fp3q5f(i,2)= double(subs(solyz.z(2,1))); % INVALID NOT VERIFIED
end
% tab1q5= horzcat((fm'),fp1q5f); % table with f,y,z
% tab2q5= horzcat((fm'),fp2q5f);
tab3q5= horzcat((fm'),fp3q5f); % INVALID NOT VERIFIED