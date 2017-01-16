function smoothed_z = UWMA(z,L)
%
% Computes the unweighted moving average (over lookback window L) of input
% data z
%

smoothed_z = zeros(size(z));
for i = 1:size(z,1)
    for j = L:size(z,2)
        smoothed_z(i,j) = mean(z(i,(j-(L-1)):j));
    end
end
