
%% Signal AN Lab 2

%% ##### Task 1 #####

% Task1.1 Task 1: The Periodogram and Averaging

% Given x[n] - 1.5x[n-1] + 0.64x[n-2] = e[n], AR(2) model
% The syntax is G = freqz(B,A,2*pi*f0);, where B=1; and A=[1 -1.5 0.64]; 
% represent the numerator and denominator polynomials, respectively. Plot the amplitude in dB, using 20*log10(abs(G));, versus the linear frequency scale f0.
%The Plot of the True Spectrum

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
% The given information to generate a periodogram
% This is Gaussian random noise variables, W.
N=1024; 
L=50;
W=randn(N,L); 
X=filter(1,[1 -1],W);
x=X(L+1:end);
f0=0:0.001:0.5;

XX=fft(x,N); 
P=XX.*conj(XX)/N;
P=P(1:N/2);
f=0:1/N:(N-1)/(2*N);
P_dB=10*log10(abs(P));
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

N=1024;
L=50;
W=randn(N,L); 
X=filter(1,[1 -1],W);
x=X(L+1:end);
k1=4;
k2=16;
M1=N/k1
M2=N/k2
xx=reshape(x,[M1,k1]);
xx=reshape(x,[M2,k2]);
XX=fft(xx,N);
PP=XX.*conj(XX)/M;
PB=mean(PP');
PB=PB(1:N/2)';



%% ##### Task 2 #####

%% Task2.1


%% Task2.2


%% Task2.3
PW=pwelch(x,window,[],2*pi*f)*2*pi;


%###Parametric Model Estimation  ###










A = [1 -2.7607 3.8106 -2.6535 0.9238];
[H,F] = freqz(1,A,[],1);
plot(F,20*log10(abs(H)))
title('Amplitude-Magnitude Response(H) vs Linear Frequency')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
rng default

x = randn(1000,1);
y = filter(1,A,x);
[Pxx,F] = pyulear(y,6,512,1); % Equal to true spectrum when order>5
hold on
plot(F,10*log10(Pxx))
legend('True Power Spectral Density','pyulear PSD Estimate')
% Comment on change of ORder:
% If order=1 The PSD is way bad then the true model order
% and if p=100 it became equal to the true order but with lot of added
% noise (fluctuations shows noise is added).




%% ##### Task 5 #####

% Task 5.1
% Compute N=8192 samples of the ARMA process:
%   x[n]-0.24*x[n-1] + 0.08*x[n-2] - 0.37*x[n-3] + 0.52*x[n-4] = e[n]+0.56*e[n-1]+0.81*e[n-2].

N=8192; L=1;
X=randn(N,L); % Gaussian random variables.
B=[1 0.56 0.81]; A=[1 -0.24 0.08 -0.37 0.52]; % Filter cofficeints (noise and o/p)
Y=filter(B,A,X); % Geenration of ARMA model using filter command
t = linspace(-pi,pi,N);
plot(t,X)
% hold on
plot(t,Y)
legend('Input Data','Filtered Data')
hold off
% f=0:0.001:0.5;
f = linspace(-pi,pi,512);
G=freqz(Y,2*pi*f);
G_dB=20*log10(abs(G));
figure(1)
plot(f,G_dB) % This is a true spectrum
title('Filter magnitude Response G vr frequency')
xlabel('Linear Frequency(Hz)')
ylabel('Magnitude Response(dB)')

%% Task 5.2####

%% Task 5.3####

N=8192; L=1;
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
%
% Task 6.1
d=randn(1024,1); e=randn(1024,1);
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
