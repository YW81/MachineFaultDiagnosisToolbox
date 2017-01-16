%--------------------------------------------------------------------------
% Name:            KF_Demo1.m
%
% Description:     Application of a Kalman filter to the task of tracking a
%                  white noise process or constant process based on noisy
%                  measurements.
%
% Author:          Brian Moore
%                  brimoor@umich.edu
%
% Date:            June 1, 2012
%--------------------------------------------------------------------------

% Process noise variance
Q = 0.001;

% Measurement noise variance
R = 0.02;

%
% Q and R estimation procedure
%
% 1. Specify a lookback window and EWMALength
EWMALength = 15;
lookbackWindow = 30;
% 2. Set use true variance flag to 'true' to reveal actual values to KF
UseTrueVariances = 'true';
%

% Number of iterations  
N = 500;

%--------------------------------------------------------------------------
% State variable initializations
%--------------------------------------------------------------------------
% True State
x = zeros(1,N);

% Apriori state estimates
x_apriori = zeros(1,N);

% Aposteriori state estimates
x_aposteriori = zeros(1,N);

% Apriori error covariance estimates
P_apriori = zeros(1,N);

% Aposteriori error covariance estimates
P_aposteriori = zeros(1,N);

% Measurements
z = zeros(1,N);

% Smoothed measurements
smoothed_z = zeros(1,N);

% Kalman Gain
K = zeros(1,N);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Initializations
%--------------------------------------------------------------------------
% True initial state
x(1) = rand;

% EWMA lambda value
lambda = 1-2/(EWMALength+1);

% First measurement
z(1) = x(1) + sqrt(Q)*randn;

% First smoothed measurement
smoothed_z(1) = z(1);

% Initial aposteriori state estimate
x_aposteriori(1:lookbackWindow) = .5;

% Initial aposteriori error covariance estimate
P_aposteriori(1:lookbackWindow) = 1;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Perform Kalman filtering
%--------------------------------------------------------------------------
for i = 2:N
    % Update true state
    x(i) = x(i-1) + sqrt(Q)*randn;

    % Update measurements
    z(i) = x(i) + sqrt(R)*randn;

    % Update smoothed measurement
    smoothed_z(i) = lambda * smoothed_z(i-1) + (1-lambda)*z(i);
    
    % Wait until lookback window is primed before applying Kalman filter
    if (i >= lookbackWindow)
        %-----------------------------------------------------------
        % Estimate noise/process covariance
        %-----------------------------------------------------------
        if strcmpi(UseTrueVariances,'true')
            Ri = R;
            Qi = Q;
        else
            [Ri Qi] = StandardCovEst(z,smoothed_z,i,lookbackWindow);
        end
        %-----------------------------------------------------------
        
        % Update Kalman filter
        [x_aposteriori(i) P_aposteriori(i)] = KalmanFilterIteration(z(i),Qi,Ri,x_aposteriori(i-1), P_aposteriori(i-1));
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------------
figure
% Plot true state
b = plot(2:N,x(2:N),'b');
hold on
% Plot estimates
r = plot(lookbackWindow:N,x_aposteriori(lookbackWindow:N),'r');
% Plot measurements
g = plot(2:N,z(2:N),'g+');
title(['Kalman Filtering: Tracking a White Noise Process with Q = ' num2str(Q) ', R = ' num2str(R)]);
legend([b r g],'True Value','KF Estimate','Measurements');
xlabel('Time')
ylabel('Value')
grid on
%--------------------------------------------------------------------------
