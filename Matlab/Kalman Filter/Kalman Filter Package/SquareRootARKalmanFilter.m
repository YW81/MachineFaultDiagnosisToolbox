function [x_kf varargout] = SquareRootARKalmanFilter(z,M,MAlen,varargin)
%--------------------------------------------------------------------------
% Syntax:       x_kf = SquareRootARKalmanFilter(z,M,MAlen,N);
%               x_kf = SquareRootARKalmanFilter(z,M,MAlen,'EWMA');
%               x_kf = SquareRootARKalmanFilter(z,M,MAlen,N,'EWMA');
%               x_kf = SquareRootARKalmanFilter(z,M,MAlen,N,'UWMA');
%               [x_kf KF] = SquareRootARKalmanFilter(z,M,MAlen,N);
%               [x_kf KF] = SquareRootARKalmanFilter(z,M,MAlen,'EWMA');
%               [x_kf KF] = SquareRootARKalmanFilter(z,M,MAlen,N,'EWMA');
%               [x_kf KF] = SquareRootARKalmanFilter(z,M,MAlen,N,'UWMA');
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
%               KF.UP_pr - Apriori error covariance estimates (U)
%               KF.DP_pr - Apriori error covariance estimates (D)
%               KF.UPa_pr - Apriori AR coefficients error covariance (U)
%               KF.DPa_pr - Apriori AR coefficients error covariance (D)
%               KF.UP_po - Aposteriori error covariance estimates (U)
%               KF.DP_po - Aposteriori error covariance estimates (D)
%               KF.UPa_po - Aposteriori AR coefficient error covariance (U)
%               KF.DPa_po - Aposteriori AR coefficient error covariance (D)
%               KF.UQ - Process covariance estimates (U)
%               KF.DQ - Process covariance estimates (D)
%               KF.UQa - AR coefficients state covariance estimates (U)
%               KF.DQa - AR coefficients state covariance estimates (D)
%               KF.R - Noise covariance estimates
%               KF.UR - Noise covariance estimates (U)
%               KF.DR - Noise covariance estimates (D)
%               KF.Ra - AR coefficients noise covariance estimates
%               KF.URa - AR coefficients noise covariance estimates (U)
%               KF.DRa - AR coefficients noise covariance estimates (D)
%
% Description:  This function performs square root Kalman auto-regressive
%               (AR) filtering on noisy input data z. The assumed system
%               model is that the noisy measurements (z) = linear
%               combination of M previous true state (x) samples + white
%               noise. This function estimates the process/noise
%               covariances of the input data at each iteration using a
%               smoothed version of z as a surrogate for the true process
%               state.
%
%               NOTE: This square root implementation of the Kalman AR
%               filter produces the same output as StandardARKalmanFilter()
%               except when the noisy measurements are poorly-conditioned,
%               in which case SquareRootARKalmanFilter() prodcues a MORE
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
    error('Input syntax error. Type ''help SquareRootARKalmanFilter'' for assistance.');
end

% AR Coefficient memory
MEM_RATE = 0.75; % 0 = No memory, 1 = Infinite memory

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
UP_apriori  = zeros(DIM*M,DIM*M,n);
DP_apriori  = zeros(DIM*M,DIM*M,n);
UPa_apriori = zeros(DIM*M,DIM*M,n);
DPa_apriori = zeros(DIM*M,DIM*M,n);

% Aposteriori error covariance estimates
UP_aposteriori  = zeros(DIM*M,DIM*M,n);
DP_aposteriori  = zeros(DIM*M,DIM*M,n);
UPa_aposteriori = zeros(DIM*M,DIM*M,n);
DPa_aposteriori = zeros(DIM*M,DIM*M,n);

% Process variance estimates
UQ  = zeros(DIM*M,DIM*M,n);
DQ  = zeros(DIM*M,DIM*M,n);
UQa = eye(DIM*M);
DQa = eye(DIM*M);

% Measurement variance estimates
R  = zeros(DIM,DIM,n);
UR  = zeros(DIM,DIM,n);
DR  = zeros(DIM,DIM,n);
Ra = zeros(DIM,DIM,n);
URa = zeros(DIM,DIM,n);
DRa = zeros(DIM,DIM,n);
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
UP_aposteriori(:,:,startIndex-1)  = eye(DIM*M);
DP_aposteriori(:,:,startIndex-1)  = eye(DIM*M);
UPa_aposteriori(:,:,startIndex-1) = eye(DIM*M);
DPa_aposteriori(:,:,startIndex-1) = eye(DIM*M);

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
Qtemp = zeros(DIM);
for i = startIndex:n
    %----------------------------------------------------------------------
    % Estimate noise/process covariances
    %----------------------------------------------------------------------
    if strcmpi(CovMethod,'Standard')
        [R(:,:,i) Qtemp] = StandardCovEst(z,smoothed_z,i,N);
    else
        [R(:,:,i) Qtemp] = EWMACovEst(z,smoothed_z,i,N,R(:,:,i-1),Qtemp);
    end
    
    % Decay AR model variance
    DQa = (1 - MEM_RATE) * DQa;
    Ra(:,:,i) = R(:,:,i) + Qtemp;
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Update Ha
    %----------------------------------------------------------------------
    for j = 1:DIM
        Ha(j,(M*(j-1)+1):(M*j)) = x_aposteriori(((j-1)*M+1):(j*M),(i-1));
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Square root Kalman filter the AR parameters
    %----------------------------------------------------------------------
    try
        % Update UD decompositions
        [URa(:,:,i) DRa(:,:,i)] = myUD(Ra(:,:,i));

        % Time update
        [a_apriori(:,i) UPa_apriori(:,:,i) DPa_apriori(:,:,i)] = thornton(a_aposteriori(:,i-1),UPa_aposteriori(:,:,i-1),DPa_aposteriori(:,:,i-1),UQa,DQa);
        
        % Decorrelate measurements
        z_ind = myUnitTriSysSol(URa(:,:,i),z(:,i),'upper');
        H_ind = myUnitTriSysSol(URa(:,:,i),Ha,'upper');

        % Measurement Update
        a_aposteriori(:,i) = a_apriori(:,i);
        UPa_aposteriori(:,:,i) = UPa_apriori(:,:,i);
        DPa_aposteriori(:,:,i) = DPa_apriori(:,:,i);
        for j = 1:DIM
            [a_aposteriori(:,i) UPa_aposteriori(:,:,i) DPa_aposteriori(:,:,i)] = bierman(z_ind(j),DRa(j,j,i),H_ind(j,:),a_aposteriori(:,i),UPa_aposteriori(:,:,i),DPa_aposteriori(:,:,i));
        end
    catch e
        %disp(e.message);
        a_aposteriori(:,i) = a_aposteriori(:,i-1);
        UPa_aposteriori(:,:,i) = UPa_aposteriori(:,:,i-1);
        DPa_aposteriori(:,:,i) = DPa_aposteriori(:,:,i-1);
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Update Ai
    %----------------------------------------------------------------------
    for j = 1:DIM
        Ai((j-1)*M+1,((j-1)*M+1):(j*M)) = a_aposteriori(((j-1)*M+1):(j*M),i)';
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Square root Kalman filter the data
    %----------------------------------------------------------------------
    try
        % Update UD decompositions
        [UQtemp DQtemp] = myUD(Qtemp);
        [UR(:,:,i) DR(:,:,i)] = myUD(R(:,:,i));
        for j = 1:DIM
            for k =1:DIM
                UQ(M*(j-1)+1,M*(k-1)+1,i) = UQtemp(j,k);
                DQ(M*(j-1)+1,M*(k-1)+1,i) = DQtemp(j,k);
            end
        end

        % Time update
        [x_apriori(:,i) UP_apriori(:,:,i) DP_apriori(:,:,i)] = thornton(x_aposteriori(:,i-1),UP_aposteriori(:,:,i-1),DP_aposteriori(:,:,i-1),UQ(:,:,i),DQ(:,:,i),Ai);

        % Decorrelate measurements
        z_ind = myUnitTriSysSol(UR(:,:,i),z(:,i),'upper');
        H_ind = myUnitTriSysSol(UR(:,:,i),H,'upper');

        % Measurement Update
        x_aposteriori(:,i) = x_apriori(:,i);
        UP_aposteriori(:,:,i) = UP_apriori(:,:,i);
        DP_aposteriori(:,:,i) = DP_apriori(:,:,i);
        for j = 1:DIM
            [x_aposteriori(:,i) UP_aposteriori(:,:,i) DP_aposteriori(:,:,i)] = bierman(z_ind(j),DR(j,j,i),H_ind(j,:),x_aposteriori(:,i),UP_aposteriori(:,:,i),DP_aposteriori(:,:,i));
        end
    catch e
        %disp(e.message);
        for j = 1:DIM
            for k = M:-1:2
               x_aposteriori((j-1)*M+k,i) = x_aposteriori((j-1)*M+(k-1),i-1);
            end
            x_aposteriori((j-1)*M+1,i) = z(j,i);
        end
        UP_aposteriori(:,:,i)  = eye(DIM*M);
        DP_aposteriori(:,:,i) = eye(DIM*M);
    end
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------

% Return user requested data
x_kf = x_aposteriori;
if (nargout == 2)
    varargout{1} = struct('x_pr',x_apriori, ...
                          'a_pr',a_apriori, ...
                          'a_po',a_aposteriori, ...
                          'UP_pr',UP_apriori, ...
                          'DP_pr',DP_apriori, ...
                          'UPa_pr',UPa_apriori, ...
                          'DPa_pr',DPa_apriori, ...
                          'UP_po',UP_aposteriori, ...
                          'DP_po',DP_aposteriori, ...
                          'UPa_po',UPa_aposteriori, ...
                          'DPa_po',DPa_aposteriori, ...
                          'UQ',UQ, ...
                          'DQ',DQ, ...
                          'UQa',UQa, ...
                          'DQa',DQa, ...
                          'R',R, ...
                          'UR',UR, ...
                          'DR',DR, ...
                          'Ra',Ra, ...
                          'URa',URa, ...
                          'DRa',DRa);
end
