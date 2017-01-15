% 带通滤波器
% x 原始信号
% fs 采样频率
% lf 低通截止
% hf 高通截止
function [xx] = mybpf(x,fs,lf,hf)
N = length(x);
len=N;
F = fft(x);
L = round(lf/fs*len);
H = round(hf/fs*len);
%dd = len-cc;
F(1:L)=0;
F(N-L+2:N)=0;
F(H:N-H+2)=0;

aa =F;
xx = ifft(aa);
end