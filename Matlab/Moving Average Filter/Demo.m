% ------------------------------------------------------------------------
%                         Demo : Moving Average Filter
%--------------------------------------------------------------------------
clear;
close all;
clc;
pos = [200   200   1000  500];
% -------------------Simulation signal generation--------------------------
fs = 1000;
t = 0:1/fs:1-1/fs;
x = cos(2*pi*10*t)+0.2*randn(size(t));
x = x';
N = length(x);
%---------------------------------Main-------------------------------------%
windowlength = 10;
output = movingmean(x,windowlength,[],[]);
%--------------------------------Result------------------------------------%
% Original waveform
figure
plot(t,x)
% De-noised signal
hold on
plot(t,output)
xlabel('Time [s]')
ylabel('Magnitude')
legend('Original signal','De-noised signal');
setfontsize(14);
set(gcf,'pos',pos);


