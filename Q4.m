%=========================================================================
% I.Y. Tanasa (1034117)       
% Aleman Zapata, R.A. (59383)
% Modeling Dynamics Project
%=========================================================================
%% QUESTION 4 DERIVE 2ND ORDER DYNAMICAL SYSTEM
% close all, clear all, clc;
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

