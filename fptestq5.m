close all, clear all, clc;
e= 10^-2; p= 0.5; q=0.05; f=0.1;% condition 1 
u=0 ; % without external excitation
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
y= 0;%double(subs((3*f)/4 - (400*f^2 - 680*f + 441)^(1/2)/80 + 21/80));
z= 0;%double(subs((400*f^2 - 680*f + 441)^(1/2)/2 - 10*f + 19/2));
testyplus(i)= double(subs(-y+fm(i)*z-(((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))*y)); % for x plus
testzplus(i)= double(subs(((((4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - y + 1)/(2*q))-z)/p)); % for x plus 
testyminus(i)= double(subs(-y+fm(i)*z-(-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))*y)); % for x minus
testzminus(i)= double(subs(((-(y + (4*q*u - 2*y + 4*q*y + y^2 + 1)^(1/2) - 1)/(2*q))-z)/p)); % for x minus
end

