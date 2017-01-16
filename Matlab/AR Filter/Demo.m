% ------------------------------------------------------------------------
%                         Demo : AR filter
%--------------------------------------------------------------------------
clear;
close all;
clc;
pos = [200   200   1000  500];
% -------------------Simulation signal generation--------------------------
fs = 1000;
t = 0:1/fs:1-1/fs;
x = cos(2*pi*10*t)+0.2*randn(size(t));
N = length(x);
%---------------------------------Main-------------------------------------%
order = 512;
a =lpc(x,order);
x1 = filter([0 -a(2:end)],1,x);
x2 = x-x1;
%--------------------------------Result------------------------------------%
% Original waveform
figure
plot(t,x)
% LPC Estimate
hold on
plot(t,x1)
xlabel('Time [s]')
ylabel('Magnitude')
legend('Original signal','LPC estimate');
setfontsize(14);
set(gcf,'pos',pos);

% % Prediction Error
% figure
% plot(t,x2)
% xlabel('Time [s]')
% ylabel('Magnitude')
% setfontsize(14);
% set(gcf,'pos',pos);

