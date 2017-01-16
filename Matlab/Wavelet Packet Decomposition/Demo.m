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
load WPD.mat;
%-----------------------Initialization of parameters----------------------%
fs = 24000;
N = length(x);
t = 0:1/fs:(N-1)/fs;
%---------------------------------Main-------------------------------------%
% Perform decomposition at level 3 of s using db1. 
[c,l] = wavedec(x,4,'db1');
%--------------------------------Result------------------------------------%
% Original waveform
figure
plot(t,x)
xlabel('Time [s]')
ylabel('Magnitude')
setfontsize(14);
set(gcf,'pos',pos);

% Decomposition result
figure
for k=1:max(size(l)-1)
    subplot(max(size(l)-1),1,k);
    if k>1
        startIndex = startIndex+l(k-1)-1;
    else
        startIndex = 1;
    end
    tempx = x(startIndex:(startIndex+l(k)-1));
    plot(0:1/fs:(length(tempx)-1)/fs,tempx);
    xlabel('Time [s]')
    ylabel('Magnitude')
    setfontsize(14);
end