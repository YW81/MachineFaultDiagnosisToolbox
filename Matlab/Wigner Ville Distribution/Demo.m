% ------------------------------------------------------------------------
%                         Demo : Wigner Ville Distribution
%--------------------------------------------------------------------------
% -------------------------------------------------------------------------
%  Run simulation first to generate suitable data x,
% then implement this demo to perform corresponding algorithm
%  Or load the existed .mat file
% -------------------------------------------------------------------------
clear;
close all;
clc;
pos = [200   200   1000  500];
% ---------------------------------Load data x-----------------------------
load WVD.mat;
%-----------------------Initialization of parameters----------------------%
fs = 24000;
x = x(1:0.1*fs);
x = x';
N = length(x);
t = 0:1/fs:(N-1)/fs;
f = 0:fs/N:fs/2;
%---------------------------------Main-------------------------------------%
[tfr, tt, ff] = wv(x, t, N);
tfr = abs(tfr);
tfr = tfr';
%--------------------------------Result------------------------------------%
% Original waveform
figure
plot(t,x);
xlabel('Time [s]')
ylabel('Magnitude')
setfontsize(14);
set(gcf,'pos',pos);

% Original Frequency Spectrum
figure
myfft(fs,x,1);
xlabel('Time [s]')
ylabel('Frequency [Hz]')
setfontsize(14);
set(gcf,'pos',pos);

% Wigner Ville Distribution
figure
specgmShow = abs(tfr(1:(N/2+1),:));
spmax = max(max(specgmShow));
spmin = min(min(specgmShow));
colormap(jet);
image(linspace(0,N,N)/fs,(0:N/2-1)/N*fs,256*(specgmShow-spmin)/(spmax-spmin));
axis('xy')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
setfontsize(14);
set(gcf,'pos',pos);
