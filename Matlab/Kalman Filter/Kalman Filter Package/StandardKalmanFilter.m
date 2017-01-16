function [x_kf varargout] = StandardKalmanFilter(z,MAlen,varargin)
%--------------------------------------------------------------------------
% Syntax:       x_kf = StandardKalmanFilter(z,MAlen,N);
%               x_kf = StandardKalmanFilter(z,MAlen,'EWMA');
%               x_kf = StandardKalmanFilter(z,MAlen,N,'EWMA');
%               x_kf = StandardKalmanFilter(z,MAlen,N,'UWMA');
%               [x_kf KF] = StandardKalmanFilter(z,MAlen,N);
%               [x_kf KF] = StandardKalmanFilter(z,MAlen,'EWMA');
%               [x_kf KF] = StandardKalmanFilter(z,MAlen,N,'EWMA');
%               [x_kf KF] = StandardKalmanFilter(z,MAlen,N,'UWMA');
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
%               KF.P_pr - Apriori error covariance estimates
%               KF.P_po - Aposteriori error covariance estimates
%               KF.K - Kalman gains
%               KF.Q - Process covariance estimates
%               KF.R - Noise covariance estimates
%
% Description:  This function performs standard Kalman filtering on noisy
%               input data z. The assumed system model is that the noisy
%               measurements (z) = the true state (x) + white noise. This
%               function estimates the process/noise covariances of the
%               input data at each iteration using a smoothed version of z
%               as a surrogate for the true process state.
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
    error('Input syntax error. Type ''help StandardKalmanFilter'' for assistance.');
end

% Kalman gain computation method
KalmanGainMethod = 'pinv'; % 'pinv' or 'UD'

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
P_apriori = zeros(DIM,DIM,n);

% Aposteriori error covariance estimates
P_aposteriori = zeros(DIM,DIM,n);

% Kalman Gain
K = zeros(DIM,DIM,n);

% Process variance estimates
Q = zeros(DIM,DIM,n);

% Measurement variance estimates
R = zeros(DIM,DIM,n);
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
P_aposteriori(:,:,startIndex-1) = eye(DIM);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Now perform Kalman filtering
%--------------------------------------------------------------------------
for i = startIndex:n
    % Update apriori estimate
    x_apriori(:,i) = x_aposteriori(:,i-1);

    % Estimate noise/process covariance
    if strcmpi(CovMethod,'Standard')
        [R(:,:,i) Q(:,:,i)] = StandardCovEst(z,smoothed_z,i,N);
    else
        [R(:,:,i) Q(:,:,i)] = EWMACovEst(z,smoothed_z,i,N,R(:,:,i-1),Q(:,:,i-1));
    end
    
    % Update apriori error covariance estimate
    P_apriori(:,:,i) = P_aposteriori(:,:,i-1) + Q(:,:,i);

    %----------------------------------------------------------------------
    % Compute Kalman Gain
    %----------------------------------------------------------------------
    if strcmpi(KalmanGainMethod,'pinv')
        % PInv method
        [inv isSingular] = myPInv(P_apriori(:,:,i) + R(:,:,i));
        K(:,:,i) = P_apriori(:,:,i) * inv;
    else
        % UD method
        [K(:,:,i) isSingular] = KalmanGainCalc(P_apriori(:,:,i),R(:,:,i));
    end
    %----------------------------------------------------------------------
    
    % Measurement Update
    if sum(sum(isnan(K(:,:,i)))) || strcmpi(isSingular,'true')
        x_aposteriori(:,i) = z(:,i);
        P_aposteriori(:,:,i) = eye(DIM);
    else
        % Update aposteriori state estimate
        x_aposteriori(:,i) = x_apriori(:,i) + K(:,:,i) * (z(:,i) - x_apriori(:,i));

        % Update aposteriori error covariance estimate
        P_aposteriori(:,:,i) = (eye(DIM) - K(:,:,i)) * P_apriori(:,:,i);
    end
end
%--------------------------------------------------------------------------

% Return user requested data
x_kf = x_aposteriori;
if (nargout == 2)
    varargout{1} = struct('x_pr',x_apriori, ...
                          'P_pr',P_apriori, ...
                          'P_po',P_aposteriori, ...
                          'K',K, ...
                          'Q',Q, ...
                          'R',R);
end
