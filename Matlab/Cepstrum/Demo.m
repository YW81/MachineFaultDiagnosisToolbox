% ------------------------------------------------------------------------
%                         Demo : Cepstrum
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
load Cepstrum.mat;
%-----------------------Initialization of parameters----------------------%
fs = 24000;
N = length(x);
t = 0:1/fs:(N-1)/fs;
%---------------------------------Main-------------------------------------%
y = real(ifft(log(abs(fft(x)))));
%--------------------------------Result------------------------------------%
% Original waveform
figure
plot(t,x)
xlabel('Time [s]')
ylabel('Magnitude')
setfontsize(14);
set(gcf,'pos',pos);

% Cepstrum
figure
plot(t,y)
xlabel('Time [s]')
ylabel('Magnitude')
setfontsize(14);
set(gcf,'pos',pos);
xlim([0 0.06])
ylim([-0.05 0.1])