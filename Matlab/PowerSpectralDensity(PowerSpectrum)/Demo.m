% ---------------------------------------------------------------
%          Demo : Power Spectral Density (Power Spectrum)
%--------------------------------------------------------------------------
%-----------------------Initialization of parameters----------------------%
clear;
close all;
clc;
pos = [200   200   1000  500];
% -------------------Simulation signal generation--------------------------
fs = 1000;
t = 0:1/fs:1-1/fs;
x = cos(2*pi*100*t)+randn(size(t));
N = length(x);
%---------------------------------Main-------------------------------------%
[pxx,f] = pwelch(x,fs);
% Normalization
pxx = pxx/(fs/2)*pi;
ff = f/pi*fs/2;
%--------------------------------Result------------------------------------%
% Original waveform
figure
plot(t,x)
xlabel('Time (s)')
ylabel('Magnitude')
setfontsize(14);
set(gcf,'pos',pos);

% Magnitude(dB)
figure
plot(ff,10*log10(pxx))
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
setfontsize(14);
set(gcf,'pos',pos);

% Magnitude
figure
plot(ff,pxx)
xlabel('Frequency [Hz]')
ylabel('Magnitude')
setfontsize(14);
set(gcf,'pos',pos);