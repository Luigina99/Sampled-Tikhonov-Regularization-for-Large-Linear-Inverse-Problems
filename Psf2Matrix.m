function [A] = Psf2Matrix(A)

n=size(A, 1);
B=zeros(n, n);

for i=1:n
    B(:, i)=A*[zeros(i-1, 1);1;zeros(n-i, 1)];
end

A=B;
