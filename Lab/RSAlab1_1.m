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

x = randn(1,2000);  % Generate x vectors
y = randn(1,2000);  % Generate y vectors

plot(x,y,'.')
title('Scatterplot x VS y, Normal Distribution')

%% 2.2 b

u = -sqrt(3)+2*sqrt(3)*rand(1,2000); % generate u vectors 
v = -sqrt(3)+2*sqrt(3)*rand(1,2000); % generate v vectors

plot(u,v,'.')                        % generate Scatterplot u VS v
title('Scatterplot u VS v, Uniform Distribution')  

plot(x,u,'.')                         % generate Scatterplot x VS u
title('Scatterplot x VS u, Uniform Distribution')  


%% 2.3
N=2000;              % 2000 points
yh=0.5;              % fixed value ŷ
dy=0.1;              % tolerance δy
xk= randn(1,N);      % x random variable with Normal distributed
yk= randn(1,N);      % y random variable with Normal distributed
k=1;                 % Realization start 
for(i=1:N)           % starting point for loop to serch
    if (yk(i)>=yh-dy)&&(yk(i)<=yh+dy)    % condition of the loop xk where ŷ − ∆y < y and k < ŷ + ∆y
        x(k)=xk(i);                      % save the value that of x then we will plot this value
        k=k+1;                           % move to next realization
    end
end

figure()
histogram(x);
title('Histogram for selected values with a fixed y')

%% 2.3 b

N=2000;           % 2000 points
uh=0.5;           % fixed value u
du=0.1;           % tolerance δu
xku= rand(1,N);   % x random variable with Normal distributed
uk= rand(1,N);    % u random variable with Normal distributed
k=1;              % Realization start 
for(i=1:N)
    if (uk(i)>=uh-du)&&(uk(i)<=uh+du)   % condition of the loop xk where uh − ∆u < u and k < uh + ∆u
        x(k)=xku(i);   % save the value that of x then we will plot this value
        k=k+1;   % move to next realization
    end
end

figure()
hist(x, 20);
title('Histogram for selected values with a fixed y')

%% 3.1 a= 0.5

 x = randn(1,2000);            % x random variable with Normal distributed on 2000 points
 y = randn(1,2000);            % y random variable with Normal distributed on 2000 points
 
 a= 0.5;                      % alpha 
 
 z = (a*x + sqrt(1-(a.^2))*y);  % z equation
 
 figure()
 plot(x,z,'.')                 % Scatterplot for x VS alpha=0.5 
 title('Scatterplot for alpha=0.5')
 
 figure()
 hist(z)
 title('Histogram for alpha=0.5') % Histogram for alpha=0.5
 
 %% 3.1 a= -0.5

 x = randn(1,2000);            % x random variable with Normal distributed on 2000 points
 y = randn(1,2000);            % y random variable with Normal distributed on 2000 points
 
 a= -0.5;                      % alpha 
 
 z = (a*x + sqrt(1-(a.^2))*y);  % z equation
  
 figure()
 plot(x,z,'.')
 title('Scatterplot for alpha=-0.5')
 
 figure()
 hist(z)
 title('Histogram for alpha=-0.5')
 
 %% 3.1 a= 0.9

 x = randn(1,2000);            % x random variable with Normal distributed on 2000 points
 y = randn(1,2000);            % y random variable with Normal distributed on 2000 points
 
 a= 0.9;                       % alpha 
 
 z = (a*x + sqrt(1-(a.^2))*y); % z equation
 
 figure()
 plot(x,z,'.')                  % Scatterplot for alpha=0.9
 title('Scatterplot for alpha=0.9')
 
 figure()
 hist(z)   % Histogram for alpha=0.9
 title('Histogram for alpha=0.9')
 %% 3.1 a= -0.9

 x = randn(1,2000);            % x random variable with Normal distributed on 2000 points
 y = randn(1,2000);            % y random variable with Normal distributed on 2000 points
 
 a= -0.9;                      % alpha 
 
 z = (a*x + sqrt(1-(a.^2))*y);  % z equation
 
 figure()
 plot(x,z,'.')   % Scatterplot for alpha=-0.9
 title('Scatterplot for alpha=-0.9')

 figure()
 hist(z)     % Histogram for alpha=-0.9
 title('Histogram for alpha=-0.9')
 
 %% 3.2 a = 0.7 and p(x|z = 0.5)
 
 x = randn(1,2000);            % x random variable with Normal distributed on 2000 points
 y = randn(1,2000);            % y random variable with Normal distributed on 2000 points
 
 a= 0.7;                       % alpha 
 
 z = (a*x + sqrt(1-(a.^2))*y);  % z equation
 
 zh = 0.5;        % distribution p(x|z = 0.5)
 dz = 0.1;        % tolerance
 
k=1;
N=2000;         % number of point

for (i=1:N)
    if (z(i) <= zh+dz && z(i) >= zh-dz)   % condition 
        xz(k) = x(i);
        k=k+1;
    end
end
figure()
hist(xz)
title('Histogram for alpha=0.7')

 %% 3.2 a = 0.7 and p(x|z = -0.5)
 
 x = randn(1,2000);            % x random variable with Normal distributed on 2000 points
 y = randn(1,2000);            % y random variable with Normal distributed on 2000 points
 
 a= 0.7;                       % alpha 
 
 z = (a*x + sqrt(1-(a.^2))*y);  % z equation
 
 zh = -0.5;        % distribution p(x|z = 0.5)
 dz = 0.1;         % tolerance
 
k=1;
N=2000;           % number of point

for (i=1:N)
    if (z(i) <= zh+dz && z(i) >= zh-dz)   % condition 
        xz(k) = x(i);
        k=k+1;
    end
end
figure()
hist(xz)
title('Histogram for alpha=0.7')

%% 3.2 a = -0.7 and p(x|z = 0.5)
 
 x = randn(1,2000);            % x random variable with Normal distributed on 2000 points
 y = randn(1,2000);            % y random variable with Normal distributed on 2000 points
 
 a= -0.7;                     % alpha 
 
 z = (a*x + sqrt(1-(a.^2))*y);  % z equation
 
 zh = 0.5;                     % distribution p(x|z = 0.5)
 dz = 0.1;                     % tolerance
 
k=1;
N=2000;        % number of point

for (i=1:N) 
    if (z(i) <= zh+dz && z(i) >= zh-dz)  % condition
        xz(k) = x(i);
        k=k+1;
    end
end

figure()
hist(xz)
title('Histogram for alpha=-0.7')
 
%% 4.1

N = 256;                        % number of row
K= 256;                         % number of culumes 
X=randn(N,K);                   % nurmal distributed matrix
Ensemble_avg = mean(X');        % Ensemble Avreg value of a sample for the matrix
Time_avg= mean(X);              % Time Averag value for the kth realization of the matrix

figure()
plot(Ensemble_avg)              % Ensemble Avreg figure
hold on
plot(Time_avg)                  % Time Averag figure
hold off

%% 4.2


N = 256;                        % number of row
K= 256;                         % number of culumes 
X=randn(N,K);                   % nurmal distributed matrix

for (i=1:5)                     % afew values for n
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
for j=1:3
    subplot(2,2,j)
    plot(X(n3(j),:),X(n4(j),:),'.')
end

%% 5.5

for i=1:256
    w = randn(1,256);
    x = filter(1,[1 -1],w);
    X(:,i) = x;
end

for n=2:256
    r(n-1) = (X(n,:)*X(n-1,:)')/256; % Average
end

plot(r)
hold on
y=(0:256);
x=y;
plot(x,y,'--')
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
    %rxt(n-1) = 0.9*((1-0.9^(n)));
end
figure()
plot(r)
hold on
%plot(rxt,'--')
hold off

