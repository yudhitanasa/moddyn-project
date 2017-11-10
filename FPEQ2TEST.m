close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; p= 0.5; q=0.05; % condition 1
f=0.1;
u=0 ; % without external excitation
% FP1 test
y= 0;
z= 0;
% FP2 test
% y= 0.580206097986845;
% z= -1.20824391947380;
% FP3 test
% y= 0.0947939020131550;
% z= 18.2082439194738;
% double(subs(-y + f*z - x*y + u))
% double(subs((x-z) / p))
testyplus= double(subs(-y+f*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))*y)) % for x plus
testzplus= double(subs(((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))-z)/p)) % for x plus 

testyminus= double(subs(-y+f*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))*y)) % for x minus
testzminus= double(subs(((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p)) % for x minus