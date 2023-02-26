
clc; clear all; close all;

%--------------- PARAMETERS ---------------  

struct.T1=1000;
struct.T2=100;
struct.TR=10;
struct.TE=(struct.TR/2);
struct.alpha=60;
struct.beta=-400:1:400; % in degrees
I=eye(3);

%---------------------------------------------

phi=[0];
p=length(phi); % number of periods p
signal=[];
TR=struct.TR;
alpha=struct.alpha;

for df=struct.beta
    if p==1 % bSSFP case
      A1= A(struct,df,TR,phi(1),alpha);
      B1=B(struct,TR);
      Mss=inv(I-A1)*B1;
      signal=[signal,Mss(1)+1i*Mss(2)];
    else % higher order 
     B_term=B(struct,TR);
     Mss=left_term(@A,struct, df,TR,phi,alpha,2)*right_term(@A,struct,df,TR,phi,alpha,p)*B_term;
     signal=[signal,Mss(1)+1i*Mss(2)];    
    end

end

plot(struct.beta,abs(signal));
grid on;
xlabel('off-resonance (degrees)');
ylabel('Amplitude');
