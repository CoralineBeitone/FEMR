clear;
clc;

T1=1000;
T2=100;
TR=10;
TE=(TR/2);
alpha=60;
beta=-400:1:400; % in degrees
I=eye(3);
M0=[0;0;1];
signal=[];

for df=beta
    R1 = Rmat(alpha, 0);
    P1 = Pmat(df, TR, TR);
    E1 = Emat(T1, T2, TR);
    A1 = P1*E1*R1;
    B1 = (I - E1) * M0;
    
    R2 = Rmat(alpha, 90);
    P2 = Pmat(df, TR, TR);
    E2 = Emat(T1, T2, TR);
    A2 = P2*E2*R2;
    B2 = (I - E2) * M0;

    ETE = Emat(T1, T2, TE);
    PTE = Pmat(df, TR, TE);


    Mneg = inv(I-A1*A2)*(A2*B1+B2);
    Mpos = A2*Mneg+B2;                                   
    MTE = PTE*ETE*Mpos + (I-ETE)*M0;

    signal = [signal MTE(1)+1j*MTE(2)];
    
end
plot(abs(signal))





function R = Rmat(alpha, phi)
    c=cosd(alpha);
    s=sind(alpha);
    x= cosd(phi);
    y=sind(phi);
    
    R=[ c*y*y+x*x, (1-c)*x*y, -s*y; 
        (1-c)*x*y, c*x*x+y*y, s*x;
        s*y, -s*x, c
        
    ];
end 

function P = Pmat(beta, TR, t)
    P=[cosd(beta*(t/TR)), sind(beta*(t/TR)), 0; -sind(beta*(t/TR)), cosd(beta*(t/TR)), 0; 0, 0 ,1];
end 

function E = Emat(T1, T2, t)
    E1=exp(-t/T1);
    E2=exp(-t/T2);
    E=[E2, 0, 0; 0, E2, 0; 0, 0 , E1];
end 
