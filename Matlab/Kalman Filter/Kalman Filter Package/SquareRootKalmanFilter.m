function [x_kf varargout] = SquareRootKalmanFilter(z,MAlen,varargin)
%--------------------------------------------------------------------------
% Syntax:       x_kf = SquareRootKalmanFilter(z,MAlen,N);
%               x_kf = SquareRootKalmanFilter(z,MAlen,'EWMA');
%               x_kf = SquareRootKalmanFilter(z,MAlen,N,'EWMA');
%               x_kf = SquareRootKalmanFilter(z,MAlen,N,'UWMA');
%               [x_kf KF] = SquareRootKalmanFilter(z,MAlen,N);
%               [x_kf KF] = SquareRootKalmanFilter(z,MAlen,'EWMA');
%               [x_kf KF] = SquareRootKalmanFilter(z,MAlen,N,'EWMA');
%               [x_kf KF] = SquareRootKalmanFilter(z,MAlen,N,'UWMA');
%
% Inputs:       z is an m x n matrix containing n samples of an
%               m-dimensional signal
%
%               MAlen is the length of the moving average to apply to z
%               (used during covariance estimation)
%
%               N is the length of the lookback window to use during
%               covariance estimation. When N is specified and nargin == 3,
%               this function applies an exponentially weighted moving
%               average (with memory N) during covariance estimation
%
%               When nargin == 3, 'EWMA' instructs this function to use
%               exponentially weighted (recursive) covariance estimation
%
%               When nargin == 4, 'EWMA' instructs this function to apply
%               an exponentially weighted moving average with (memory N) to
%               z during covariance estimation. Alternatively, 'UWMA'
%               applies an unweighted moving average of length N to z
%
% Outputs:      x_kf is an m x n matrix containing the Kalman filter
%               "true state" estimate derived from noisy measurements z.
%
%               KF is a struct containing the following fields:
%
%               KF.x_pr - Apriori state estimates
%               KF.UP_pr - Apriori error covariance estimates (U)
%               KF.DP_pr - Apriori error covariance estimates (D)
%               KF.UP_po - Aposteriori error covariance estimates (U)
%               KF.DP_po - Aposteriori error covariance estimates (D)
%               KF.Q - Process covariance estimates
%               KF.UQ - Process covariance estimates (U)
%               KF.DQ - Process covariance estimates (D)
%               KF.R - Process covariance estimates
%               KF.UR - Process covariance estimates (U)
%               KF.DR - Noise covariance estimates (D)
%
% Description:  This function performs square root Kalman filtering on
%               noisy input data z. The assumed system model is that the
%               noisy measurements (z) = the true state (x) + white noise.
%               This function estimates the process/noise covariances of
%               the input data at each iteration using a smoothed version
%               of z as a surrogate for the true process state.
%
%               NOTE: This square root implementation of the Kalman filter
%               produces the same output as StandardKalmanFilter(), except
%               when the noisy measurements are poorly-conditioned, in
%               which case SquareRootKalmanFilter() prodcues a MORE
%               NUMERICALLY STABLE output (at the cost of greater
%               computational complexity.)
%
%               NOTE: "Square root" refers to the fact that, instead of
%               computing the error covariance P at each iteration, this
%               function computes U and D such that [U D] = myUD(P). This
%               strategy is known to be more numerically stable.
%
% Author:       Brian Moore
%               brimoor@umich.edu
%
% Date:         July 12, 2012
%--------------------------------------------------------------------------

% Parse user inputs
if (nargin == 3)
    if ~ischar(varargin{1})
        N = varargin{1};
        MAType = 'EWMA';
        CovMethod = 'Standard';
    else
        N = MAlen;
        MAType = 'EWMA';
        CovMethod = 'EWMA';
    end
elseif (nargin == 4)
    N = varargin{1};
    MAType = varargin{2};
    CovMethod = 'Standard';
else
    error('Input syntax error. Type ''help SquareRootKalmanFilter'' for assistance.');
end

%--------------------------------------------------------------------------
% State variable initializations
%--------------------------------------------------------------------------
% Data dimension
DIM = size(z,1);

% Number of iterations
n = size(z,2);

% Apriori state estimates
x_apriori = zeros(DIM,n);

% Aposteriori state estimates
x_aposteriori = zeros(DIM,n);

% Apriori error covariance estimates
UP_apriori = zeros(DIM,DIM,n);
DP_apriori = zeros(DIM,DIM,n);

% Aposteriori error covariance estimates
UP_aposteriori = zeros(DIM,DIM,n);
DP_aposteriori = zeros(DIM,DIM,n);

% Process variance estimates
Q = zeros(DIM,DIM,n);
UQ = zeros(DIM,DIM,n);
DQ = zeros(DIM,DIM,n);

% Measurement variance estimates
R = zeros(DIM,DIM,n);
UR = zeros(DIM,DIM,n);
DR = zeros(DIM,DIM,n);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% More initializations
%--------------------------------------------------------------------------
% MA smoothed measurements
eval(['smoothed_z = ' MAType '(z,MAlen);']);

% Simulation start index
startIndex = N + MAlen - 1;

% Initial aposteriori state estimate
x_aposteriori(:,startIndex-1) = smoothed_z(:,startIndex-1); %#ok

% Initial aposteriori error covariance estimate
UP_aposteriori(:,:,startIndex-1) = eye(DIM);
DP_aposteriori(:,:,startIndex-1) = eye(DIM);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Now perform square root Kalman filtering
%--------------------------------------------------------------------------
for i = startIndex:n
    % Estimate noise/process covariance
    if strcmpi(CovMethod,'Standard')
        [R(:,:,i) Q(:,:,i)] = StandardCovEst(z,smoothed_z,i,N);
    else
        [R(:,:,i) Q(:,:,i)] = EWMACovEst(z,smoothed_z,i,N,R(:,:,i-1),Q(:,:,i-1));
    end
    try
        % Update UD decompositions
        [UQ(:,:,i) DQ(:,:,i)] = myUD(Q(:,:,i));
        [UR(:,:,i) DR(:,:,i)] = myUD(R(:,:,i));
        
        % Time update
        [x_apriori(:,i) UP_apriori(:,:,i) DP_apriori(:,:,i)] = thornton(x_aposteriori(:,i-1),UP_aposteriori(:,:,i-1),DP_aposteriori(:,:,i-1),UQ(:,:,i),DQ(:,:,i));
        
        % Decorrelate measurements
        z_ind = myUnitTriSysSol(UR(:,:,i),z(:,i),'upper');
        H_ind = myUnitTriSysSol(UR(:,:,i),eye(DIM),'upper');

        % Measurement Update
        x_aposteriori(:,i) = x_apriori(:,i);
        UP_aposteriori(:,:,i) = UP_apriori(:,:,i);
        DP_aposteriori(:,:,i) = DP_apriori(:,:,i);
        for j = 1:DIM
            [x_aposteriori(:,i) UP_aposteriori(:,:,i) DP_aposteriori(:,:,i)] = bierman(z_ind(j),DR(j,j,i),H_ind(j,:),x_aposteriori(:,i),UP_aposteriori(:,:,i),DP_aposteriori(:,:,i));
        end
    catch e
        %disp(e.message);
        x_aposteriori(:,i) = z(:,i);
        UP_aposteriori(:,:,i)  = eye(DIM);
        DP_aposteriori(:,:,i) = eye(DIM);
    end 
end
%--------------------------------------------------------------------------

% Return user requested data
x_kf = x_aposteriori;
if (nargout == 2)
    varargout{1} = struct('x_pr',x_apriori, ...
                          'UP_pr',UP_apriori, ...
                          'DP_pr',DP_apriori, ...
                          'UP_po',UP_aposteriori, ...
                          'DP_po',DP_aposteriori, ...
                          'Q',Q, ...
                          'UQ',UQ, ...
                          'DQ',DQ, ...
                          'R',R, ...
                          'UR',UR, ...
                          'DR',DR);
end
