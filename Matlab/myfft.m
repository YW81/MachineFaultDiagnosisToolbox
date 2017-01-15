%  信号作FFT变换，得到信号的频谱
%  本程序不适用于宽带信号

function [ff,amp] = myfft(fs,x,draw,linewidth,style)

if nargin  < 5
    style = '-';
end

if nargin  == 2
    draw = 0;
end

if nargin  == 3
    linewidth = 1;
end

NN = length(x);

% ff = (0:floor(NN/2))*fs/NN;
% amp = abs(fft(x)/NN*2);
% amp = amp(1:(floor(NN/2)+1));

ff = linspace(0,fs,NN+1);
ff = ff(1:NN);
amp = abs(fft(x)/NN*2);

amp = amp(1:round(NN/2));
ff = ff(1:round(NN/2));

if draw == 1
    plot(ff,amp,style,'linewidth',linewidth);
end

end