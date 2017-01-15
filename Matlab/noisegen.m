function [Y,NOISE] = noisegen(X,SNR)
% noisegen add white Gaussian noise to a signal.
% [Y, NOISE] = NOISEGEN(X,SNR) adds white Gaussian NOISE to X.  The SNR is in dB.
NOISE=randn(size(X));
NOISE=NOISE-mean(NOISE);
signal_power = 1/length(X)*sum(X.*X);
noise_variance = signal_power / ( 10^(SNR/10) );
NOISE=sqrt(noise_variance)/std(NOISE)*NOISE;
Y=X+NOISE;