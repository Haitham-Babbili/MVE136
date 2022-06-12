
%% Signal AN Lab 2

%% ##### Task 1 #####

% Task1.1

B=1; 
A=[1,-1.5,0.64]; %Filter cofficients, from the AR(2) model
f0 =[0:0.001:0.5] ;
G=freqz(B,A,2*pi*f0);
G_dB=20*log10(abs(G));
figure(1)
plot(f0,G_dB) % This is the plot of the true spectrum
title('Amplitude-Magnitude Response(G) vs Linear Frequency')
xlabel('Frequency (Hz)')
ylabel('Magnitude(dB)')
grid on


%% Task 1.2
% Generate a periodogram This is Gaussian noise , W.
B=1; 
A=[1,-1.5,0.64];

N=1024; 
L=50;
W=randn(1,N+L); 
x=filter(1,[1 -1.5 0.64],W);
x=x(L+1:end);
f0=0:0.001:0.5;

X=fft(x,N); 
P=X.*conj(X)/N;
P=P(1:N/2);
P_dB=10*log10(abs(P));
f=0:1/N:(N-1)/(2*N);
G=freqz(B,A,2*pi*f0);
G_dB=20*log10(abs(G));
figure(2)
plot(f,P_dB)
title('Periodogram (dB) vs Frequency (Hz)')
xlabel('Frequency (Hz)')
ylabel('Periodogram(dB)')
grid on

hold on
plot(f0,G_dB)
hold off

%% Task 1.3

B=1; 
A=[1,-1.5,0.64];
N=1024;
L=50;
f0=0:1/N:(N-1)/(2*N);
W=randn(1,N+L); 
x=filter(1,[1 -1.5 0.64],W);
x=x(L+1:end);
k1=4;
k2=16;
M1=N/k1;
M2=N/k2;

xx1=reshape(x,M1,k1);
XX1=fft(xx1,N);
PP1=XX1.*conj(XX1)/M1;
PB1=mean(PP1');
PB1=PB1(1:N/2)';
PB1_dB=10*log10(abs(PB1));

xx2=reshape(x,M2,k2);
XX2=fft(xx2,N);
PP=XX2.*conj(XX2)/M2;
PB=mean(PP');
PB=PB(1:N/2)';
PB_dB=10*log10(abs(PB));

figure(3)
% plot(XX1)
% hold on 
% plot(XX2)
% hold off
figure(3)
plot(f0,PB1_dB)
hold on 

plot(f0,PB_dB)
hold off

%% ##### Task 2 #####

%% Task2.1
wintool;

%% Task2.2


%% Task2.3

PW = pwelch(x,window,[],2*pi*f);





%% ##### Task 3 #####

% Task 3.1

%x=;
%M=;
%NFFT=;



[PBT,fgrid]= btmethod(x,M,NFFT);

%% Task 3.2
%x=;
M=24;
NFFT=1024;

[PBT,fgrid]= btmethod(x,M,NFFT);
plot(PBT)

%% ##### Task 4 #####

A = [1 -2.7607 3.8106 -2.6535 0.9238];
%f_ar=(NFFT/2)+1;
[H,f_ar] = freqz(1,A,[],1);
plot(f_ar,20*log10(abs(H)))
title('Amplitude-Magnitude Response(H) vs Linear Frequency')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
rng default

x = randn(1000,1);
y = filter(1,A,x);
[PAR,f_ar] = pyulear(y,6,512,1); % Equal to true spectrum when order>5
hold on
plot(f_ar,10*log10(PAR))
legend('True Power Spectral Density','pyulear PSD Estimate')




%% ##### Task 5 #####

% Task 5.1
%x[n]-0.24*x[n-1] + 0.08*x[n-2] - 0.37*x[n-3] + 0.52*x[n-4] = e[n]+0.56*e[n-1]+0.81*e[n-2].

N=8192; 
L=1;
X=randn(N,L); % Gaussian random variables.
B=[1 0.56 0.81];  % Filter cofficeints (noise and o/p)
A=[1 -0.24 0.08 -0.37 0.52]; % Filter cofficeints (noise and o/p)
Y=filter(B,A,X); % create of ARMA model using filter command
% f1 = linspace(-pi,pi,N);
% plot(f1,X)
% hold on
% plot(f1,Y)
% legend('Input Data','Filtered Data')
% hold off
% f=0:0.001:0.5;
f = linspace(-pi,pi,512);
G=freqz(Y,2*pi*f);
G_dB=20*log10(abs(G));
figure(10)
plot(f,G_dB) % This is a true spectrum
title('Filter magnitude Response G vr frequency')
xlabel('Linear Frequency(Hz)')
ylabel('Magnitude Response(dB)')

%% Task 5.2

%% Task 5.3

N=8192; 
L=1;
X=randn(N,L); % Gaussian random variables.
B=[1]; A=[1 -0.24 0.08 -0.37 0.52]; % Filter cofficeints (noise and o/p)
Y=filter(B,A,X); % Geenration of AR model using filter command
% True spectrum
f=0:0.001:0.5;
G=freqz(Y,2*pi*f);
G_dB=20*log10(abs(G));
plot(G_dB)
hold on
[Pxx,F] = pyulear(Y,25,1024,1); 
plot(F,10*log10(Pxx))
legend('True Power Spectral Density','pyulear PSD Estimate')







%% ##### Task 6 #####

% x[n]=d[n]+w[n];
%   d[n]-0.13*d[n-1] + 0.9*d[n-2]= ed[n] (AR Process)
%   w[n]=ew[n]-0.8ew*d[n-1] + 0.2*ew[n-2] (MA Process)
% AR filter cofficeients
B_ar=[1]; 
A_ar=[1 -0.13 0.9]; 
B_ma=[1 -0.8 0.2]; 
A_ma=[1];

% Task 6.1
d=randn(1,1024); 
e=randn(1,1024);
Y_ar=filter(B_ar,A_ar,d); 
Y_ma=filter(B_ma,A_ma,e); 
f=0:0.001:0.5;
G=freqz(Y_ar,2*pi*f);
G_dB_AR=20*log10(abs(G));
plot(G_dB_AR)
hold on    
G=freqz(Y_ma,2*pi*f);
G_dB_MA=20*log10(abs(G));
plot(G_dB_MA)
title('Filter magnitude Response G vr frequency')
xlabel('Linear frequency')
ylabel('Magnitude response')
legend('Desired Signal response','Disturbance response')
% Comments: The graphs shows that the desired signal can select a range of
% frequencies which is a band pass characteristics whereas noise act as
% high pass filter because its response value started increasing with the
% increase in frequency which proves that it act is a high pass filter.
