function [result] = left_term(A,param, beta,t,phi,alpha,p)

I=eye(3);
product=I;

for j=1:p
    product=A(param,beta,t,phi(j),alpha)*product;
end

result=(I-product);

if det(result)==0
    print('Matrix non invertible');
else
    result=inv(result);
end


end

