function [K isSingular] = KalmanGainCalc(P,R,varargin)
%
% Computes the Kalman Gain K from the formula K = P*H'*(H*P*H'+R)^(-1)
% WITHOUT explicitly inverting H*P*H'+R by computing its unit Cholesky
% Decomposition and then applying backsubstitutions on the resulting
% triangular factors.
%

if nargin ~= 3
    try
        [U D] = myUD(P+R);
        isSingular = 'false';
        X1 = myUnitTriSysSol(U,P','upper');
        X2 = X1;
        for i = 1:size(X2,1)
            for j = 1:size(X2,2)
                X2(i,j) = X2(i,j) / D(i,i);
            end
        end
        X3 = myUnitTriSysSol(U',X2,'lower');
        K = X3';
    catch %#ok
        isSingular = 'true';
        K = zeros(size(P,1),size(R,2));
    end
else
    H = varargin{1};
    try
        [U D] = myUD(H*P*H'+R);
        isSingular = 'false';
        X1 = myUnitTriSysSol(U,H*P','upper');
        X2 = X1;
        for i = 1:size(X2,1)
            for j = 1:size(X2,2)
                X2(i,j) = X2(i,j) / D(i,i);
            end
        end
        X3 = myUnitTriSysSol(U',X2,'lower');
        K = X3';
    catch %#ok
        isSingular = 'true';
        K = zeros(size(P,1),size(R,2));
    end
end
