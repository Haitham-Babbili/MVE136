
%% ##### Task 1 #####
% Task 1.1

f0=0:0.001:0.5;
A = [1 -1.5 0.64];
B = 1;

G = freqz(B,A,2*pi*f0);
y = 20*log10(abs(G));

figure(1)
plot(f0,y)
xlabel('Frequency (Hz)')
ylabel('Amplitude dB')
title('True Spectrum')

%% Task 1.2
% Generate a periodogram This is Gaussian noise , W.
N = 1024;
L = 50;

b = randn(1,N+L);
A = [1 -1.5 0.64];
f=0:1/N:(N-1)/(2*N); 

x = filter(1,A,b);

x=x(L+1:end);

X=fft(x,N);
P=X.*conj(X)/N;
P=P(1:N/2);
PdB=10*log10(abs(P));

f0=0:0.001:0.5;
G = freqz(1,A,2*pi*f0);
y = 20*log10(abs(G));

figure(2)
plot(f,PdB)
 hold on
 plot(f0,y)
 hold off
 title('Periodogram (dB) vs Frequency (Hz)')
xlabel('Frequency (Hz)')
ylabel('Periodogram(dB)')
 %% Task 1.3
 
 N = 1024;
 L = 50;
 K1 = 4;
 M1 = N/K1;
 f=0:1/N:(N-1)/(2*N);
 b = randn(1,N+L);
 A = [1 -1.5 0.64];
 
 x = filter(1,A,b);
 x=x(L+1:end);
 
 
 xx1=reshape(x,M1,K1);
 XX1=fft(xx1,N);
 PP1=XX1.*conj(XX1)/M1;
 PB1=mean(PP1');
 PB1=PB1(1:N/2)';
 
 plot(f,10*log10(abs(PB1)))
 hold on
 
 K2 = 16;
 M2 = N/K2;
 
 xx2=reshape(x,M2,K2);
 XX2=fft(xx2,N);
 PP2=XX2.*conj(XX2)/M2;
 PB2=mean(PP2');
 PB2=PB2(1:N/2)';
 
 figure(3)
 plot(f,10*log10(abs(PB2)))
 hold on
 
 f0=0:0.001:0.5;
 G = freqz(1,A,2*pi*f0);
 y = 20*log10(abs(G));
 

 plot(f0,y)
 hold off
 
 title('Periodogram Average')
 legend('Bartlett, K=4','Bartlett, K=16','True Spectrum')
 xlabel('Frequency (Hz)')
ylabel('Periodogram(dB)')
 
 %% ##### Task 2 #####
 
 %% Task 2.1
 wintool;

 %% Task 2.2

 %% Task 2.3
 
 N = 1024;
 L = 50;
 K1 = 4;
 M1 = N/K1;
 f=0:1/N:(N-1)/(2*N);
 b = randn(1,N+L);
 A = [1 -1.5 0.64];
 
x = filter(1,A,b);
x=x(L+1:end);
w= window(@hamming,M1);
figure(4)
PW = pwelch(x,w,[],2*pi*f)*2*pi;
plot(f,10*log10(abs(PW)))
hold on

 % Pure estimation
 K2 = 16;
 M2 = N/K2;
 
 xx2=reshape(x,M2,K2);
 XX2=fft(xx2,N);
 PP2=XX2.*conj(XX2)/M2;
 PB2=mean(PP2');
 PB2=PB2(1:N/2)';
 %figure(77)
 plot(f,10*log10(abs(PB2)))
 hold off
 
 title('Periodogram Average Vs Pure estimation of Periodogram Average')
 legend('Welch','Bartlett, K=16')
 xlabel('Frequency (Hz)')
 ylabel('Periodogram(dB)')
 
 %% ##### Task 3 #####
 
 % Task 3.1
 
N = 1024;
L = 50;
b = randn(1,N+L);
x = filter(1,[1 -1.5 0.64], b);
x = x(L+1:end);
M1 = 30;
%M1 = 10;
NFFT = 1024;
 
[PBT, fgrid] = btmethod(x, M1, NFFT);
figure(5)
plot(fgrid,PBT)
hold on
%  for M1=10:10:40
%     NFFT = 1024;
%     [PBT, fgrid] = btmethod(x, M1, NFFT);
%     figure(5)
%     plot(fgrid,PBT)
%     hold on
%  end


% Task 3.2 Bartlett

K = 16;
f=0:1/N:(N-1)/(2*N);
b = randn(1,N+L);
A = [1 -1.5 0.64];
 
x = filter(1,A,b);
x=x(L+1:end);

M = N/K;
xx = reshape(x,M,K);
XX = fft(xx,N);
PP = XX.*conj(XX)/M;
PB = mean(PP');
PB = PB(1:N/2);
 
 plot(f,10*log10(abs(PB)));

 hold off
 title('Estemated Spectrum Vs Bartlett ')
 legend('Blackman-Tukey','Bartlett')
 xlabel('Frequency (Hz)')
 ylabel('Periodogram(dB)')
 
%% ##### Task 4 #####
% Task 4.1
N = 1024;
L = 50;
b = randn(1,N+L);
x = filter(1,[1 -1.5 0.64], b);
x = x(L+1:end);

order = 2;

[PAR,f_ar] = pyulear(x,order,N);
f_ar=f_ar/(2*pi);
PAR=PAR.*pi;

figure(6)
plot(f_ar,10*log10(abs(PAR)));
hold on
% True spectrum
f0=0:0.001:0.5;
A = [1 -1.5 0.64];
B = 1;

G = freqz(B,A,2*pi*f0);
y = 20*log10(abs(G));

plot(f0,y)
hold off

legend('Parametric AR','True Spectrum');
title('Parametric AR Modeling')
xlabel('Frequency (Hz)')
ylabel('Periodogram(dB)')
%% Task 4.2

N = 1024;
L = 50;
b = randn(1,N+L);
x = filter(1,[1 -1.5 0.64], b);
x = x(L+1:end);

order = 1;

[PAR,f_ar] = pyulear(x,order,N);
f_ar=f_ar/(2*pi);
PAR=PAR.*pi;

figure(7)
plot(f_ar,10*log10(abs(PAR)));
hold on

N = 1024;
L = 50;
b = randn(1,N+L);
x = filter(1,[1 -1.5 0.64], b);
x = x(L+1:end);

order = 100;

[PAR,f_ar] = pyulear(x,order,N);
f_ar=f_ar/(2*pi);
PAR=PAR.*pi;


plot(f_ar,10*log10(abs(PAR)));
hold on

% True spectrum
f0=0:0.001:0.5;
A = [1 -1.5 0.64];
B = 1;

G = freqz(B,A,2*pi*f0);
y = 20*log10(abs(G));

plot(f0,y)
hold off
title('Parametric AR Modeling')
legend('Parametric AR, p=1','Parametric AR, p=100','True Spectrum')
xlabel('Frequency (Hz)')
ylabel('Periodogram(dB)')
%% ##### Task 5 #####

% 5.1

N = 8192;
f0=0:0.001:0.5;

A = [1 -0.24 0.08 -0.37 0.52];
B = [1 0.56 0.81];

G = freqz(B,A,2*pi*f0);
Y = 20*log10(abs(G));

figure(8)
plot(f0, Y)
title('True Spectrum, ARMA')
xlabel('Frequency (Hz)')
ylabel('Periodogram(dB)')
%% 5.2
NFFT = 2048;

%% 5.3

N = 8192;
A = [1 -0.24 0.08 -0.37 0.52];
B = [1 0.56 0.81];

random = randn(1,N);
x = filter(B,A,random);

NFFT = 2048;

for p=10:10:50
    
    [PAR,f_ar] = pyulear(x,p,NFFT);
    f_ar=f_ar/(2*pi);
    PAR=PAR.*pi;
    figure(10)
    plot(f_ar,10*log10(abs(PAR)));
    hold on
end

% True Spectrum
f0=0:0.001:0.5;
G = freqz(B,A,2*pi*f0);
Y = 20*log10(abs(G));

plot(f0, Y,'Linewidth',2)

hold off
title('AR Modeling')
legend('Order = 10','Order = 20','Order = 30','Order = 40','Order = 50','True')
xlabel('Frequency (Hz)')
ylabel('Periodogram(dB)')

%% ##### Task 6 #####
% 6.1

f0=0:0.001:0.5;

A_d = [1 -0.13 0.9];
B_d = 1;

A_w = 1;
B_w = [1 -0.8 0.2];

G_d = freqz(B_d,A_d,2*pi*f0);
Y_d = 10*log10(abs(G_d));

G_w = freqz(B_w,A_w,2*pi*f0);
Y_w = 10*log10(abs(G_w));

figure(11)
plot(f0, Y_d)
hold on
plot(f0,Y_w)
hold off

title('Filter magnitude Response G vr frequency')
legend('Desired Signal response','Disturbance response')
legend('Desired[n]','Disturbance[n]')
xlabel('Linear frequency')
ylabel('Magnitude response')


%% 6.A.3

figure(12)

% p=1
r_d1 = 5.2879;
r_w1 = 1.68;
r_x1 = r_d1 + r_w1;
R_X1 = toeplitz(r_x1);
h_1 = R_X1\r_d1;
[H1,W1] = freqz(h_1);
subplot(2,2,1)
plot(W1/pi,20*log10(abs(H1)));
title('p=1')

% p=2
r_d2 = [5.2879 ;0.3618];
r_w2 = [1.68; -0.96];
r_x2 = r_d2 + r_w2;
R_X2 = toeplitz(r_x2);
h_2 = R_X2\r_d2;
[H2,W2] = freqz(h_2);
subplot(2,2,2)
plot(W2/pi,20*log10(abs(H2)));
title('p=2')

% p=5
r_d5(1) = 5.2879;
r_d5(2) = 0.3618;
r_w5 = [1.68; -0.96; 0.2; 0; 0]; 
for k = 3:5
    r_d5(k) = 0.13.*r_d5(k-1) - 0.9.*r_d5(k-2);   
end
r_d5 = r_d5';
r_x5 = r_d5 + r_w5;
R_X5 = toeplitz(r_x5);
h_5 = R_X5\r_d5;
[H5,W5] = freqz(h_5);
subplot(2,2,3)
plot(W5/pi,20*log10(abs(H5)));
title('p=5')

% p=10
r_d10(1) = 5.2879;
r_d10(2) = 0.3618;
r_w10 = [1.68; -0.96; 0.2; 0; 0; 0; 0; 0; 0; 0]; 
for k = 3:10
    r_d10(k) = 0.13.*r_d10(k-1) - 0.9.*r_d10(k-2);   
end
r_d10 = r_d10';
r_x10 = r_d10 + r_w10;
R_X10 = toeplitz(r_x10);
h_10 = R_X10\r_d10;
[H10,W10] = freqz(h_10);
subplot(2,2,4)
plot(W10/pi,20*log10(abs(H10)));
title('p=10')



%% 6.B

% h = []
% p = [1 2 5 10];
% 
% for k = 1:4
%     sum(k) = 0;
%     for i = 1:p(k)
%         sum(k) = 

% h_1
%h = [length(h_1) lenthg(h_2) length(h_5) length(h_10)]
%p = length(h_1);
%r_d = zeros(p,1);
%r_w = zeros(p,1);
% if p == 1
%     r_d(1) = 5.2879;
%     r_w(1) = 1.68;
% elseif p == 2
%     r_d = [5.2879;0.3618];
%     r_w = [1.68;-0.96];
% elseif p >= 3
%     r_d(1) = 5.2879;
%     r_d(2) = 0.3618;
%     for k = 3:p
%         r_d(k) = 0.13*r_d(k-1)-0.9*r_d(k-2);
%     end
%     r_w = [1.68;-0.96;0.2;zeros(p-3,1)];
% end
%     r_x = r_d + r_w;
%     R_x = toeplitz(r_x);
    
    r_d_x1 = R_X1*h_1;
mse_1 = 5.2879 - sum(r_d_x1.*h_1)

% h_2
    r_d_x2 = R_X2*h_2;
mse_2 = 5.2879 - sum(r_d_x2.*h_2)

% h_5
    r_d_x5 = R_X5*h_5;
mse_5 = 5.2879 - sum(r_d_x5.*h_5)

% h_10
    r_d_x10 = R_X10*h_10;
mse_10 = 5.2879 - sum(r_d_x10.*h_10)

%% 6.B.5

N = 200;

ed = randn(1,N);
ew = randn(1,N);

A_d = [1 -0.13 0.9];
B_d = 1;

A_w = 1;
B_w = [1 -0.8 0.2];

d_n = filter(B_d,A_d,ed);
w_n = filter(B_w,A_w,ew);

x_n = d_n + w_n;

plot(d_n(150:end))
hold on
plot(x_n(150:end))

 % p=1
d_pn = conv(x_n,h_1);
d_pn = d_pn(1:N);
e_p1 = d_n-d_pn;
MSE_emp1 = norm(e_p1)^2/N

% % p=2
% d_pn2 = conv(x_n,h_2);
% d_pn2 = d_pn2(1:N);
% e_p2 = d_n2-d_pn2;
% MSE_emp2 = norm(e_p2)^2/N
% 
% % p=5
% d_pn5 = conv(x_n,h_5);
% d_pn5 = d_pn5(1:N);
% e_p5 = d_n5-d_pn5;
% MSE_emp5 = norm(e_p5)^2/N
% 
% % p=10
% d_pn10 = conv(x_n,h_10);
% d_pn10 = d_pn10(1:N);
% e_p10 = d_n10-d_pn10;
% MSE_emp10 = norm(e_p10)^2/N
% 
% 











