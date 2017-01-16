% ---------------------------------------------------------------
%                               Demo : HNR-gram
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
load HNR-gram.mat;
%-----------------------Initialization of parameters----------------------%
nlevel = 5;
global fs;
fs = 24000;
global N ;
N = length(x);
t = 0:1/fs:(N-1)/fs;
%------------------------------------Main---------------------------------%
[Bw,fc] = HNR_gram(x,nlevel,fs,1);
%-----------------------------------Result -------------------------------%
% Filter generation
if fc-Bw/2==0
    b = fir1(512,2*fc/fs,'low');
else if fc+Bw/2==fs
        b = fir1(512,2*fc/fs,'high');
    else
        b = fir1(512,2*[fc-Bw/2 fc+Bw/2]/fs);
    end
end

tempx = filtfilt(b,1,x);
hilx = abs(hilbert(tempx))-mean(abs(hilbert(tempx)));

% 最优滤波频带滤波后的时域波形 Waveform after filtration
figure
plot(t,hilx);
ylabel('Amplitude');
xlabel('Time [s]');
setfontsize(14);
set(gcf,'pos',pos);

% 最优滤波频带滤波后的包络谱  Frequency spectrum after filtration
figure
myfft(fs,hilx,1);
ylabel('Amplitude');
xlabel('Frequency [Hz]');
setfontsize(14);
set(gcf,'pos',pos);
xlim([0 300])
