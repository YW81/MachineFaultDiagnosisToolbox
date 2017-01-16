function smoothed_z = EWMA(z,L)
%
% Computes the exponentially weighted moving average (with memory L) of
% input data z
%

lambda = 1-2/(L+1);

smoothed_z = zeros(size(z));
for i = 1:size(z,1)
    smoothed_z(i,1) = z(i,1);
    for j = 2:size(z,2)
        smoothed_z(i,j) = lambda * smoothed_z(i,j-1) + (1-lambda) * z(i,j);
    end
end
