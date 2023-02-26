function [M] = B(param,t)

T1=param.T1;
TR=param.TR;
T2=param.T2;
M0=[0;0;1];

E1=exp(-t/T1);
E2=exp(-t/T2);

I=eye(3);

E=[E2, 0, 0; 0, E2, 0; 0, 0 , E1];

M=(I-E)*M0;

end

