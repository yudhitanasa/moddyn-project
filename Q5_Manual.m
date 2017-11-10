close all, clear all, clc;
syms e f p q u x y z  % redefine
% e= 10^-2; p= 0.5; q=0.05; % condition 1
% u=0 ; % without external excitation
eq22= -y+f*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))*y == 0; % 2nd order system y dot, with z=x
eq23= ((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))-z)/p == 0; % 2nd order system z dot
solyz= solve([eq22, eq23 ], [y,z]) % two set of fixed point (y,z)
% e= 10^-2; p= 0.5; q=0.05; % condition 1
% u=0 ; % without external excitation
eq22m= -y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))*y == 0; % 2nd order system y dot, with z=x
eq23m= ((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p == 0; % 2nd order system z dot
solyzm= solve([eq22m, eq23m ], [y,z]) % two set of fixed point (y,z)
e= 10^-2; p= 0.5; q=0.05; % condition 1
u=0 ; % without external excitation
jacobian([eq22 eq23],[y,z])
return



















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
