
clc
clear all
%[Noyono,fs]=audioread('vlc-record-2016-09-17-23h00m28s-Noyono_demo.mp3-.mp3');
data=load('S001.txt');                              % loading 1-channel data

Sampling_time_period=2;                             % sampling time period

sampling_frequency=500;                             % sampling frequency

[n,m]=size(data);                                   % obtain size of data

t=(1:n)*Sampling_time_period;                       % generates time vector

h=figure

plot(t,data,'p');
title('EEG Data')
grid on                                             % plot with grid
h1=figure
plot(t,data(:,1), 'b-')                             % plot without grid


%-----------------------------------------------------------------------%
%POWER SPECTRUM USING FFT
%-----------------------------------------------------------------------%


y=fft(data);                                        % fast fourier transform
power_spectrum=abs(y).^2;                           % power spectrum using fft
frequency=(1:n)*sampling_frequency/n;               % frequency vector
h2=figure
plot(frequency,20*log(power_spectrum),'b')          % plotting frequency
title('Power spectrum using fft')


%------------------------------------------------------%
%SPECTROGRAM of channel 1
%------------------------------------------------------%


[S1,F,T] = spectrogram(data(:,1),chebwin(128,100),0,sampling_frequency);

%--------------------------------------------------------------------------------%

%w=chebwin(L,r)returns the column vector w containing length L chebyshev window ;
%r is sidelobe attenuation ;
%gives 128 point chebyshev window. r = 100 dB here
 
%--------------------------------------------------------------------------------%

S1=abs(S1);
h4=figure
mesh(T,F,S1);
xlabel('Time (sec)','FontSize',14);
ylabel('Frequency (Hz)','FontSize',14);
zlabel('S1','FontSize',14);

%-------------------------------------------------------------------------%
%FREQUENCY AND SINGLE SIDED AMPLITUDE SPECTRUMS OF APLHA BETA GAMA DELTA
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%DELTA BAND-PASS FILTER (0-4)
%-------------------------------------------------------------------------%

sampling_frequency = 500;               % Sampling Frequency

Fpass = 0;                              % Passband Frequency

Fstop = 4;                              % Stopband Frequency

Dpass = 0.057501127785;                 % Passband Ripple

Dstop = 0.0001;                         % Stopband Attenuation

density  = 20;                          % Density Factor

% Calculating the order from the parameters using FIRPMORD.

[n, fo, ao, w] = firpmord([Fpass, Fstop]/(sampling_frequency/2), [1 0], [Dpass, Dstop]); 

% fo is the frequency ao is the amplitude response vector

% Calculating the coefficients using the FIRPM function.

b1 = firpm(n, fo, ao, w, {density});

Hd1 = dfilt.dffir(b1);                                 %discrete time direct form fir filter
x1=filter(Hd1,data);
h9=figure
plot(t,x1,'k')
title('waveform for DELTA band')

%================================================================%
%FREQUENCY SPECTRUM OF DELTA
%================================================================%

L=10;
sampling_frequency=500;
N_point_fft = 2^nextpow2(L);                           % Next power of 2 from length of x1
y1 = fft(x1,N_point_fft)/L;                            % for N point fft
f = sampling_frequency/2*linspace(0,1,N_point_fft/2);  % N_point_fft/2 points b/w 0 and 1

%================================================================%
%SINGLE SIDED AMPLITUDE SPECTRUM
%================================================================%

h10=figure
plot(f,2*abs(y1(1:N_point_fft/2))) 
title('Single-Sided Amplitude Spectrum of DELTA x1(t)')
xlabel('Frequency (Hz)')
ylabel('|Y1(f)|')



%-----------------------------------------------------------------%
%THETA- BAND PASS FILTER (4-7)
%-----------------------------------------------------------------%


sampling_frequency = 500; % Sampling Frequency

Fstop1 = 3.5;             % First Stopband Frequency

Fpass1 = 4;               % First Passband Frequency

Fpass2 = 7;               % Second Passband Frequency

Fstop2 = 7.5;             % Second Stopband Frequency

Dstop1 = 0.0001;           % First Stopband Attenuation

Dpass  = 0.057501127785;  % Passband Ripple

Dstop2 = 0.0001;          % Second Stopband Attenuation

density   = 20;           % Density Factor

% Calculating the order from the parameters using FIRPMORD.

[n, fo, ao, w] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(sampling_frequency/2), [0 1 0], [Dstop1 Dpass Dstop2]);

% fo is the frequency ao is the amplitude response vector

% Calculating the coefficients using the FIRPM function.

b2 = firpm(n, fo, ao, w, {density});
Hd2 = dfilt.dffir(b2);                                  %discrete time direct form fir filter
x2=filter(Hd2,data);
h11=figure
plot(t,x2,'r')
title('waveform for THETA band')


%===================================================%
%FREQUENCY SPECTRUM OF THETA 
%===================================================%


L=10;
sampling_frequency=500;
N_point_fft = 2^nextpow2(L);                            % Next power of 2 from length of x2
y2 = fft(x2,N_point_fft)/L;                             % for N point fft
f = sampling_frequency/2*linspace(0,1,N_point_fft/2);   % N_point_fft/2 points b/w 0 and 1


%===================================================%
%SINGLE SIDED AMPLITUDE SPECTRUM 
%===================================================%


h12=figure
plot(f,2*abs(y2(1:N_point_fft/2))) 
title('Single-Sided Amplitude Spectrum of THETA x2(t)')
xlabel('Frequency (Hz)')
ylabel('|Y2(f)|')


%----------------------------------------------------%
%ALPHA BAND PASS FILTER (8-12)
%----------------------------------------------------%

sampling_frequency = 500; % Sampling Frequency

Fstop1 = 7.5;             % First Stopband Frequency

Fpass1 = 8;               % First Passband Frequency

Fpass2 = 12;              % Second Passband Frequency

Fstop2 = 12.5;            % Second Stopband Frequency

Dstop1 = 0.0001;          % First Stopband Attenuation

Dpass  = 0.057501127785;  % Passband Ripple

Dstop2 = 0.0001;          % Second Stopband Attenuation

density   = 20;           % Density Factor

% Calculating the order from the parameters using FIRPMORD.

[n, fo, ao, w] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(sampling_frequency/2), [0 1 0], [Dstop1 Dpass Dstop2]);

% fo is the frequency ao is the amplitude response vector

% Calculating the coefficients using the FIRPM function.

b3  = firpm(n, fo, ao, w, {density});
Hd3 = dfilt.dffir(b3);                                  %discrete time direct form fir filter
x3=filter(Hd3,data);
h13=figure
plot(t,x3,'g')
title('waveform for ALPHA band')

%=========================================================%
%FREQUENCY SPECTRUM OF ALPHA BAND
%=========================================================%

L=10;
sampling_frequency=500;
N_point_fft = 2^nextpow2(L);                            % Next power of 2 from length of x3
y3 = fft(x3,N_point_fft)/L;                             % for N point fft
f = sampling_frequency/2*linspace(0,1,N_point_fft/2);   % N_point_fft/2 points b/w 0 and 1


%==========================================================%
%SINGLE SIDED AMPLITUDE SPECTRUM
%==========================================================%


h14=figure
plot(f,2*abs(y3(1:N_point_fft/2))) 
title('Single-Sided Amplitude Spectrum of ALPHA x3(t)')
xlabel('Frequency (Hz)')
ylabel('|Y3(f)|')


%------------------------------------------------%
%BETA  BAND PASS FILTER (12-30)
%------------------------------------------------%


sampling_frequency = 500; % Sampling Frequency

Fstop1 = 11.5;            % First Stopband Frequency

Fpass1 = 12;              % First Passband Frequency

Fpass2 = 30;              % Second Passband Frequency

Fstop2 = 30.5;            % Second Stopband Frequency

Dstop1 = 0.0001;          % First Stopband Attenuation

Dpass  = 0.057501127785;  % Passband Ripple

Dstop2 = 0.0001;          % Second Stopband Attenuation

density   = 20;           % Density Factor

% Calculating the order from the parameters using FIRPMORD.

[n, fo, ao, w] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(sampling_frequency/2), [0 1 0], [Dstop1 Dpass Dstop2]);

% fo is the frequency ao is the amplitude response vector

% Calculatng the coefficients using the FIRPM function

b4   = firpm(n, fo, ao, w, {density});

Hd4 = dfilt.dffir(b4);                                  %discrete time direct form fir filter
x4=filter(Hd4,data);
h15=figure
plot(t,x4,'y')
title('waveform for BETA band')


%==========================================================%
%FREQUENCY SPECTRUM OF BETA
%==========================================================%


L=10;
sampling_frequency=500;
N_point_fft = 2^nextpow2(L);                            % Next power of 2 from length of x4
y4 = fft(x4,N_point_fft)/L;                             % for N point fft
f = sampling_frequency/2*linspace(0,1,N_point_fft/2);   % N_point_fft/2 points b/w 0 and 1


%==========================================================%
%SINGLE SIDED AMPLITUDE SPECTRUM
%==========================================================%


h16=figure
plot(f,2*abs(y4(1:N_point_fft/2))) 
title('Single-Sided Amplitude Spectrum of BETA x4(t)')
xlabel('Frequency (Hz)')
ylabel('|Y4(f)|')


%============================================================%
%AMPLITUDE SPECTRUMS TOGETHER
%============================================================%

h17=figure
plot(f,2*abs(y1(1:N_point_fft/2)),'k',f,2*abs(y2(1:N_point_fft/2)),'r',f,2*abs(y3(1:N_point_fft/2)),'g',f,2*abs(y4(1:N_point_fft/2)),'b')
title('Single-Sided Amplitude Spectrums')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend('DELTA','THETA','APLHA','BETA');










