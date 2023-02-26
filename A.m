function [M] = A(param,beta,t,phi,alpha)


T1=param.T1;
TR=param.TR;
T2=param.T2;

c=cosd(alpha);
s=sind(alpha);
x= cosd(phi);
y=sind(phi);

R=[ c*y*y+x*x, (1-c)*x*y, -s*y; 
    (1-c)*x*y, c*x*x+y*y, s*x;
    s*y, -s*x, c
    
];


P=[cosd(beta*(t/TR)), sind(beta*(t/TR)), 0; -sind(beta*(t/TR)), cosd(beta*(t/TR)), 0; 0, 0 ,1];

E1=exp(-t/T1);
E2=exp(-t/T2);

E=[E2, 0, 0; 0, E2, 0; 0, 0 , E1];

M=P*E*R;


end

