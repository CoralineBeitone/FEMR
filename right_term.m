function [result] = right_term(A, param, beta,t,phi,alpha,p)

I=eye(3);

sum=0*I;

product=I;

assert(length(phi)==p,'Check dimension of vector phi');

if p==2
    product= A(param,beta,t,phi(p),alpha)*product;
    sum=product;

else
    for k=2:p-1
        product=I;
        for j=k:p
            product=A(param,beta,t,phi(j),alpha)*product;
        end
        sum=sum+product;
    end
end

result=sum+I;
