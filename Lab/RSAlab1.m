% 2.1
figure()
N = randn(1,2000); % generates 2000 points that are Normal distributed
hist(N,10)
title('Histogram for Normal Distribution')


figure()
U = rand(1,2000); % generates 2000 points that are Uniform distributed
hist(U,10)
title('Histogram for Uniform Distribution')

%% 2.2 a

x = randn(1,2000);
y = randn(1,2000);

plot(x,y,'.')
title('Scatterplot, Normal Distribution')

%% 2.2 b

u = -sqrt(3)+2*sqrt(3)*rand(1,2000);
v = -sqrt(3)+2*sqrt(3)*rand(1,2000);

plot(u,v,'.');
title('Scatterplot, Uniform Distribution')
hold on

figure()
plot(x,u,'.')
title('Scatterplot X versus U, Uniform Distribution') % Stander Normal (0,1)

%% 2.3
N=2000;
yh=0.5;
dy=0.1;
xk= randn(1,N);
yk= randn(1,N);
k=1;
for(i=1:N)
    if (yk(i)>=yh-dy)&&(yk(i)<=yh+dy)
        x(k)=xk(i);
        k=k+1;
    end
end


figure()
histogram(x);
title('Histogram for selected values with a fixed y')

 figure()
 histogram(yk);
 title('Histogram for random values with for xk')
% no they are not the same becuse the xk is all probabiltiy of x with
% respect of all y while xk is the probability of x withe respecte of
% smalle interval to y 
 
%figure()
%histogram(yk(i),1);
%title('Histogram for selected values with a fixed sdsd')

%% 2.3 b

N=2000;
uh=0.5;
du=0.1;
xku= rand(1,N);
uk= rand(1,N);
k=1;
for(i=1:N)
    if (uk(i)>=uh-du)&&(uk(i)<=uh+du)
        x(k)=xku(i);
        k=k+1;
    end
end

figure()
hist(x, 20);
title('Histogram for selected values with a fixed y')

%% 3.1 a= 0.5

 x = randn(1,2000);
 y = randn(1,2000);
 
 a= 0.5;
 
 z = (a*x + sqrt(1-(a.^2))*y);
 
 figure()
 plot(x,z,'.')
 title('Scatterplot for alpha=0.5')
 
 figure()
 hist(z)
 title('Histogram for alpha=0.5')
 
 %% 3.1 a= -0.5

 x = randn(1,2000);
 y = randn(1,2000);
 
 a= -0.5;
 
 z = (a*x + sqrt(1-(a.^2))*y);
 
 figure()
 plot(x,z,'.')
 title('Scatterplot for alpha=-0.5')
 
 figure()
 hist(z)
 title('Histogram for alpha=-0.5')
 
 %% 3.1 a= 0.9

 x = randn(1,2000);
 y = randn(1,2000);
 
 a= 0.9;
 
 z = (a*x + sqrt(1-(a.^2))*y);
 
 figure()
 plot(x,z,'.')
 title('Scatterplot for alpha=0.9')
 
 figure()
 hist(z)
 title('Histogram for alpha=0.9')
 %% 3.1 a= -0.9

 x = randn(1,2000);
 y = randn(1,2000);
 
 a= -0.9;
 
 z = (a*x + sqrt(1-(a.^2))*y);
 
 figure()
 plot(x,z,'.')
 title('Scatterplot for alpha=-0.9')

 figure()
 hist(z)
 title('Histogram for alpha=-0.9')
 % Here We can relize clearly at alpha = -9 the mean value
 % is 0
 
 %% 3.2 a = 0.7 and zh = 0.5
 
 x = randn(1,2000);
 y = randn(1,2000);
 
 a= 0.7;
 
 z = (a*x + sqrt(1-(a.^2))*y);
 
 zh = 0.5;
 dz = 0.1;
 
k=1;
N=2000;

for (i=1:N)
    if (z(i) <= zh+dz && z(i) >= zh-dz)
        xz(k) = x(i);
        k=k+1;
    end
end
figure()
hist(xz,30)
title('Histogram for alpha=0.7 and p(xz)=0.5')

 %% 3.2 a = 0.7 and zh = -0.5
 
 x = randn(1,2000);
 y = randn(1,2000);
 
 a= 0.7;
 
 z = (a*x + sqrt(1-(a.^2))*y);
 
 zh = -0.5;
 dz = 0.1;
 
k=1;
N=2000;

for (i=1:N)
    if (z(i) <= zh+dz && z(i) >= zh-dz)
        xz(k) = x(i);
        k=k+1;
    end
end
figure()
hist(xz,30)
title('Histogram for alpha=0.7  and p(xz)=-0.5')

%% 3.2 a = -0.7 and zh = 0.5
 
 x = randn(1,2000);
 y = randn(1,2000);
 
 a= -0.7;
 
 z = (a*x + sqrt(1-(a.^2))*y);
 
 zh = 0.5;
 dz = 0.1;
 
k=1;
N=2000;

for (i=1:N)
    if (z(i) <= zh+dz && z(i) >= zh-dz)
        xz(k) = x(i);
        k=k+1;
    end
end

figure()
hist(xz,30)
title('Histogram for alpha=-0.7  and p(xz)=0.5')
%the center is shifted by - but they have same varines   
%% 4.1

N = 256;
K= 256;
X=randn(N,K);
AVG = mean(X');
figure()
% plot(AVG)
% title('Ensemble Average')
% 
AVG2= mean(X);
% figure()
% plot(AVG2);
% title('Time Average')
plot(AVG)
hold on
plot(AVG2)
hold off
% they are ergodic in the mean because of the limet of the 

%% 4.2


N = 256;
K= 256;
X=randn(N,K);

for (i=1:2)
    n1 = randi([1,N]);
    n2 = randi([1,N]);
    figure()
    plot(X(n1,:),X(n2,:),'.')
end

%% 5.4

for i=1:256
    w = randn(1,256);
    x = filter(1,[1 -1],w);
    X(:,i) = x;
end
figure()
plot(X);
title('All Realizations')

n1 = [10, 50, 100, 200];
n2 = [9, 49, 99, 199];
figure()
for j=1:4
    subplot(2,2,j)
    plot(X(n1(j),:),X(n2(j),:),'.')
end

figure()
n3 = [50, 100, 200];
n4 = [40, 90, 190];
for i=1:3
    subplot(2,2,i)
    plot(X(n3(i),:),X(n4(i),:),'.')
end

%% 5.5

for i=1:256
    w = randn(1,256);
    x = filter(1,[1 -1],w);
    X(:,i) = x;
end

for n=2:256
    r(n-1) = (X(n,:)*X(n-1,:)')/256;
end

plot(r)
hold on
y=(0:256);
x=y;
plot(x,y,'--') %  g?r denna ber?kningen
hold off
title('Sample Ensemble Auto-correlation')
legend('Actual value','Theoretical value')
%% 6.3
for (i = 1:256)
    w = randn(1,256);
    x=filter(1,[1 -0.9],w);
    X(:,i)= x;
end

figure()
plot(X)
title('All Realizations')

n1 = [10, 50, 100, 200];
n2 = [9, 49, 99, 199];

n3 = [50, 100, 200];
n4 = [40, 90, 190];

figure()
for j=1:4
    subplot(2,2,j)
    plot(X(n1(j),:),X(n2(j),:),'.')
end

figure()
for i=1:3
    subplot(2,2,i)
    plot(X(n3(i),:),X(n4(i),:),'.')
end

%% 6.4

for i=1:256
    w = randn(1,256);
    x = filter(1,[1 -0.9],w);
    X(:,i) = x;
end

l=1;

for n=2:256
    r(n-1) = (X(n,:)*X(n-1,:)')/256;
    rxt(n-1) = 0.9*((1-0.9^(n)))/0.18;
end
figure()
plot(r)
hold on
plot(rxt,'--')
hold off

