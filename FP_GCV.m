function [x_opt] = FP_GCV(A, b, N)

A=Psf2Matrix(A);
%questo comando serve solo a trasformare la matrice A da Psf a matrice
%double

n=size(A, 1);
d=size(A, 1)/10;

%creazione matrici di campionamento
for j=1:d
    
    W(:, (j-1)*10+1:j*10)=[zeros(10*(j-1), 10); eye(10); zeros(10*(d-j), 10)];
    
end



%permutating the sampling matrices
W=RandM(W, d);


%initializing some data%
Lambda_sum=0;


A_sum=0;
A_sum_2=0;
l=10; %size of the sampling matrix%
x=ones(n, 1);

%array to plot 
X=zeros(1, N);


for i=0:N-1
    
    iterazione=i

    j=mod(i, 10)+1;
    
    if(mod(i, d)==0)       %Ciclicamente, ogni 10 it circa permutiamo la matrice W%
        W=RandM(W, d);
       end
 
    
    %calcolo elementi nella sommatoria
    W_k=W(:,(j-1)*10+1:j*10);
    
    A_k=transpose(W_k)*A;
    b_k=transpose(W_k)*b;
    A_sum=A_sum+transpose(A_k)*A_k;
    A_sum_2=A_sum_2+transpose(A_k)*transpose(W_k);
    
    %definisco le funzioni C_k e x_k in lambda
    %che mi serviranno per calcolare lambda in Dp_k
    
    C_k=@(y)((y+Lambda_sum)*eye(n)+A_sum)^(-1)*A_sum_2;
    
    x_k=@(y) C_k(y)*b;
    
    GCV_k=@(y) l*norm(A_k*x_k(y)-b_k)^2/(l-trace(A_k*C_k(y)*W_k))^2;
    
    [lambda_k, ~]=fminbnd(GCV_k, 0 , 10^(-2), optimset('Display','iter'));
    
    
     %if(i<3)
      % xp=-1:0.01:1;
       % yp=Plot(xp, GCV_k);
        %plot(xp, yp);
        %hold on
    %end
    
    Lambda_sum=Lambda_sum+lambda_k;
    
    %ora che ho trovato lambda_k, posso calcolare la k-esima iterazione x_k
    x=x_k(0);
    
   
    
end


%verificare che x_k tende a x(lambda_sum*10/N)
x_opt=x;

