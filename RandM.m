function [W] = RandM(W, d)

v=randperm(d);

for j=1:d
    
    B(:, (j-1)*10+1:j*10)= W(:,(v(j)-1)*10+1:v(j)*10);
    
end
W=B;
end