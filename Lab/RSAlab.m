% 2.1
U = rand(1,2000);
histogram(U,10)

%% 2.2 a

x = randn(1,2000);
y = randn(1,2000);

plot(x,y,'.')

%% 2.2 b

u = rand(-sqrt(3),sqrt(3));
v = rand(-sqrt(3),sqrt(3));

scatterplot(u,v,'.')
% g?r resten!!

%% 2.3
N= 2000;
yhat = 0.5;
deltay = 0.1;
yk = randn(1,N);
xk = randn(1,N);
%z=zeros(2000,1);


k = 1;

for i=1:N
    if(yk<=(yhat+deltay) & yk>=(yhat-deltay))
        x(i) = xk(k);
        k=k+1;
        i=i+1;
    end
end

histogram(x)
%hold on
R = randn(1,2000);
%plot(R,'.')
%histogram(R) %r?d
%hold off

%% 2.3 b

N= 2000;
uhat = 0.5;
deltau = 0.1;
uk = rand(1,N);
xk = rand(1,N);

k = 1;

for (i=1:N)
    if(uk<=(uhat+deltau) & uk>=(uhat-deltau))
        z(i) = xk(k);
        k=k+1;
        i=i+1;
    end
end

 histogram(z)
% hold on
U = rand(1,2000);
%plot(U,'.')
%histogram(U) %r?d
%hold off
        


