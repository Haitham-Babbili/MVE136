% 2.1
data = randn(1,2000);
histogram(data,10)
U = rand(1,2000);
histogram(U,10)

%%2.2 a
x= randn(1,2000);
y=randn(1,2000);

plot(x,y,'.');


%2.2.b
%v= randn(1,sqrt (3));u = rand(-sqrt(3),sqrt(3));
u = rand(-sqrt(3),sqrt(3));
v = rand(-sqrt(3),sqrt(3));

plot(u,v,'.')


%%2.3
yhat= 0.5;
dely=0.5;

