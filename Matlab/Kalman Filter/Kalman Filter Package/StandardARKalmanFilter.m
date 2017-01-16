function [x_kf varargout] = StandardARKalmanFilter(z,M,MAlen,varargin)
%--------------------------------------------------------------------------
% Syntax:       x_kf = StandardARKalmanFilter(z,M,MAlen,N);
%               x_kf = StandardARKalmanFilter(z,M,MAlen,'EWMA');
%               x_kf = StandardARKalmanFilter(z,M,MAlen,N,'EWMA');
%               x_kf = StandardARKalmanFilter(z,M,MAlen,N,'UWMA');
%               [x_kf KF] = StandardARKalmanFilter(z,M,MAlen,N);
%               [x_kf KF] = StandardARKalmanFilter(z,M,MAlen,'EWMA');
%               [x_kf KF] = StandardARKalmanFilter(z,M,MAlen,N,'EWMA');
%               [x_kf KF] = StandardARKalmanFilter(z,M,MAlen,N,'UWMA');
%
% Inputs:       z is an m x n matrix containing n samples of an
%               m-dimensional signal
%
%               M is the desired model order of the auto-regressive (AR)
%               model used to model the state transition.
%
%               MAlen is the length of the moving average to apply to z
%               (used during covariance estimation)
%
%               N is the length of the lookback window to use during
%               covariance estimation. When N is specified and nargin == 4,
%               this function applies an exponentially weighted moving
%               average (with memory N) during covariance estimation
%
%               When nargin == 4, 'EWMA' instructs this function to use
%               exponentially weighted (recursive) covariance estimation
%
%               When nargin == 5, 'EWMA' instructs this function to apply
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
%               KF.a_pr - Apriori AR coefficients
%               KF.a_po - Aposteriori AR coefficients
%               KF.P_pr - Apriori error covariance estimates
%               KF.Pa_pr - Apriori AR coefficients error covariance
%               KF.P_po - Aposteriori error covariance estimates
%               KF.Pa_po - Aposteriori AR coefficients error covariance
%               KF.K - Kalman gains
%               KF.Ka - Kalman gains for the AR coefficients
%               KF.Q - Process covariance estimates
%               KF.Qa - AR coefficients state covariance estimates
%               KF.R - Noise covariance estimates
%               KF.Ra - AR coefficients noise covariance estimates
%
% Description:  This function performs Kalman auto-regressive (AR)
%               filtering on noisy input data z. The assumed system model
%               is that the noisy measurements (z) = linear combination of 
%               M previous true state (x) samples + white noise. This
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
if (nargin == 4)
    if ~ischar(varargin{1})
        N = varargin{1};
        MAType = 'EWMA';
        CovMethod = 'Standard';
    else
        N = MAlen;
        MAType = 'EWMA';
        CovMethod = 'EWMA';
    end
elseif (nargin == 5)
    N = varargin{1};
    MAType = varargin{2};
    CovMethod = 'Standard';
else
    error('Input syntax error. Type ''help StandardARKalmanFilter'' for assistance.');
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
x_apriori = zeros(DIM*M,n);
a_apriori = zeros(DIM*M,n);

% Aposteriori state estimates
x_aposteriori = zeros(DIM*M,n);
a_aposteriori = zeros(DIM*M,n);

% Apriori error covariance estimates
P_apriori  = zeros(DIM*M,DIM*M,n);
Pa_apriori = zeros(DIM*M,DIM*M,n);

% Aposteriori error covariance estimates
P_aposteriori  = zeros(DIM*M,DIM*M,n);
Pa_aposteriori = zeros(DIM*M,DIM*M,n);

% Kalman Gain
K  = zeros(DIM*M,DIM,n);
Ka = zeros(DIM*M,DIM,n);

% Process variance estimates
Q  = zeros(DIM*M,DIM*M,n);
Qa = eye(DIM*M);

% Measurement variance estimates
R  = zeros(DIM,DIM,n);
Ra = zeros(DIM,DIM,n);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% More initializations
%--------------------------------------------------------------------------
% MA smoothed measurements
eval(['smoothed_z = ' MAType '(z,MAlen);']);

% Simulation start index
startIndex = N+MAlen-1;

% Initial aposteriori state estimates
for i = 1:DIM
    for j = 1:M
        for k = 1:M
            x_aposteriori((i-1)*M+j,startIndex-k) = smoothed_z(i,startIndex-(j-1)-k);
        end
    end
    a_aposteriori((i-1)*M+1,startIndex-1) = 1;
end

% Initial aposteriori error covariance estimate
P_aposteriori(:,:,startIndex-1)  = eye(DIM*M);
Pa_aposteriori(:,:,startIndex-1) = eye(DIM*M);

% Model constants
Ai = zeros(DIM*M);
for j = 1:DIM
    for k = 2:M
        Ai((j-1)*M+k,(j-1)*M+k-1) = 1;
    end
end
H = zeros(DIM,DIM*M);
for j = 1:DIM
    H(j,(j-1)*M+1) = 1;
end
Ha = zeros(DIM,DIM*M);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Perform Kalman filtering
%--------------------------------------------------------------------------
MEM_RATE = 0.75; % 0 = No memory, 1 = Infinite memory
Qtemp = zeros(DIM);
for i = startIndex:n
    % Estimate noise/process covariance
    if strcmpi(CovMethod,'Standard')
        [R(:,:,i) Qtemp] = StandardCovEst(z,smoothed_z,i,N);
    else
        [R(:,:,i) Qtemp] = EWMACovEst(z,smoothed_z,i,N,R(:,:,i-1),Qtemp);
    end
    for j = 1:DIM
        for k =1:DIM
            Q(M*(j-1)+1,M*(k-1)+1,i) = Qtemp(j,k);
        end
    end
    
    % Decay AR model variance
    Qa = (1 - MEM_RATE) * Qa;
    
    %----------------------------------------------------------------------------------------------
    % Kalman filter the AR parameters
    %----------------------------------------------------------------------------------------------
    % Update apriori estimate
    a_apriori(:,i) = a_aposteriori(:,i-1);
    % Update measurement covariance estimate
    Ra(:,:,i) = R(:,:,i) + Qtemp;
    % Update apriori error covariance estimate
    Pa_apriori(:,:,i) = Pa_aposteriori(:,:,i-1) + Qa;
    for j = 1:DIM
        Ha(j,(M*(j-1)+1):(M*j)) = x_aposteriori(((j-1)*M+1):(j*M),(i-1));
    end
    
    %----------------------------------------------------------------------
    % Compute Kalman gain
    %----------------------------------------------------------------------
    if strcmpi(KalmanGainMethod,'pinv')
        % PInv method
        [inva isSingulara] = myPInv(Ha * Pa_apriori(:,:,i) * Ha' + Ra(:,:,i));
        Ka(:,:,i) = Pa_apriori(:,:,i) * Ha' * inva;
    else
        % UD method
        [Ka(:,:,i) isSingulara] = KalmanGainCalc(Pa_apriori(:,:,i),Ra(:,:,i),Ha);
    end
    %----------------------------------------------------------------------
    
    if sum(sum(isnan(Ka(:,:,i)))) || strcmpi(isSingulara,'true')
        % Keep previous AR coefficients when inversion fails
        a_aposteriori(:,i) = a_apriori(:,i);
        Pa_aposteriori(:,:,i) = Pa_apriori(:,:,i);
    else
        % Update aposteriori state estimate
        a_aposteriori(:,i) = a_apriori(:,i) + Ka(:,:,i) * (z(:,i) - Ha * a_apriori(:,i));
        % Update aposteriori error covariance estimate
        Pa_aposteriori(:,:,i) = (eye(DIM*M) - Ka(:,:,i) * Ha) * Pa_apriori(:,:,i);
    end
    %----------------------------------------------------------------------------------------------
    
    % Update apriori estimate
    for j = 1:DIM
        Ai((j-1)*M+1,((j-1)*M+1):(j*M)) = a_aposteriori(((j-1)*M+1):(j*M),i)';
    end
    x_apriori(:,i) = Ai * x_aposteriori(:,i-1);
    
    % Update apriori error covariance estimate
    P_apriori(:,:,i) = Ai * P_aposteriori(:,:,i-1) * Ai' + Q(:,:,i);

    %----------------------------------------------------------------------
    % Update Kalman gain
    %----------------------------------------------------------------------
    if strcmpi(KalmanGainMethod,'pinv')
        % PInv method
        [inv isSingular] = myPInv(H * P_apriori(:,:,i) * H' + R(:,:,i));
        K(:,:,i) = P_apriori(:,:,i) * H' * inv;
    else
        % UD Method
        [K(:,:,i) isSingular] = KalmanGainCalc(P_apriori(:,:,i),R(:,:,i),H);
    end
    %----------------------------------------------------------------------
    
    if sum(sum(isnan(K(:,:,i)))) || strcmpi(isSingular,'true')
        K(:,:,i) = H';
    end
    % Update aposteriori state estimate
    x_aposteriori(:,i) = x_apriori(:,i) + K(:,:,i) * (z(:,i) - H * x_apriori(:,i));
    % Update aposteriori error covariance estimate
    P_aposteriori(:,:,i) = (eye(DIM*M) - K(:,:,i) * H) * P_apriori(:,:,i);
end
%--------------------------------------------------------------------------

% Return user requested data
x_kf = x_aposteriori;
if (nargout == 2)
    varargout{1} = struct('x_pr',x_apriori, ...
                          'a_pr',a_apriori, ...
                          'a_po',a_aposteriori, ...
                          'P_pr',P_apriori, ...
                          'Pa_pr',Pa_apriori, ...
                          'P_po',P_aposteriori, ...
                          'Pa_po',Pa_aposteriori, ...
                          'K',K, ...
                          'Ka',Ka, ...
                          'Q',Q, ...
                          'Qa',Qa, ...
                          'R',R, ...
                          'Ra',Ra);
end
