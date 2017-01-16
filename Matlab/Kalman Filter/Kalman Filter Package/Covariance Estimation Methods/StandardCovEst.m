function [R Q] = StandardCovEst(z,smoothed_z,i,N)
%
% Estimates the process covariance (Q) and noise covariance (R) of the data
% (z) at the current index (i) over a lookback window (N) using the
% smoothed data (smoothed_z) as a surrogate for the true process states.
%

DIM = size(z,1);

R = zeros(DIM);
Q = zeros(DIM);

if (DIM == 1)
    z = z(i-(N-1):i);
    smoothed_z = smoothed_z(i-(N-1):i);
    diffs = z - smoothed_z;

    mean1 = mean(smoothed_z);
    mean2 = mean(diffs);
    
    Q = sum((smoothed_z-mean1).^2) / (N-1);
    R = sum((diffs-mean2).^2) / (N-1);
else
    z = z(:,i-(N-1):i);
    smoothed_z = smoothed_z(:,i-(N-1):i);
    diffs = z - smoothed_z;

    for i = 1:DIM
        mean1i = mean(smoothed_z(i,:));
        mean2i = mean(diffs(i,:));
        for j = 1:DIM
            mean1j = mean(smoothed_z(j,:));
            Q(i,j) = sum((smoothed_z(i,:)-mean1i).*(smoothed_z(j,:)-mean1j))/(N-1);

            mean2j = mean(diffs(j,:));
            R(i,j) = sum((diffs(i,:)-mean2i).*(diffs(j,:)-mean2j))/(N-1);
        end
    end
end
