function [R Q] = EWMACovEst(z,smoothed_z,i,N,Rold,Qold)
%
% Updates the process covariance (Q) and noise covariance (R) estimates of
% the data (z) at the current index (i) with memory length (N) using an
% exponentially weighted averaging of the current covariances and the
% previous process covariance (Qold) and noise covariance (Rold) estimates.
% The smoothed version (smoothed_z) of the data is used as a surrogate for
% the true process states.
%

if (i == (2*N-1))
    [R Q] = StandardCovEst(z,smoothed_z,i,N);
else
    DIM = size(z,1);
    lambda = 1-2/(N+1);

    if (DIM == 1)
        Qnew = smoothed_z(i) - mean(smoothed_z((i-(N-1)):i));
        Rnew = z(i) - smoothed_z(i);
    else
        Qnew = smoothed_z(:,i) - mean(smoothed_z(:,(i-(N-1)):i),2);
        Rnew = z(:,i) - smoothed_z(:,i);
    end

    Q = (1 - lambda) * (Qnew * Qnew') + lambda * Qold;
    R = (1 - lambda) * (Rnew * Rnew') + lambda * Rold;
end
