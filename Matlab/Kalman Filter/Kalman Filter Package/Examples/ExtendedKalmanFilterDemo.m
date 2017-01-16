%--------------------------------------------------------------------------
% Name:            ExtendedKalmanFilterDemo.m
%
% Description:     Application of an Extended Kalman filter to the task of
%                  tracking a sine wave based on noisy measurments.
%
% Author:          Brian Moore
%                  brimoor@umich.edu
%
% Date:            June 1, 2012
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Model Constants
%--------------------------------------------------------------------------
% State-error Jacobian
W = [ 1 0 ;
      0 0 ];

% Measurement-state Jacobian
H = [ 1 0 ;
      0 1 ];

% Measurement-error Jacobian
V = [ 1 0 ;
      0 1 ];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% State variable initializations
%--------------------------------------------------------------------------
% Number of iterations  
N = 1000;

% True State
x = zeros(2,N);

% Apriori state estimates
x_apriori = zeros(2,N);

% Aposteriori state estimates
x_aposteriori = zeros(2,N);

% Apriori error covariance estimates
P_apriori = zeros(2,2,N);

% Aposteriori error covariance estimates
P_aposteriori = zeros(2,2,N);

% Measurements
z = zeros(2,N);

% Kalman Gain
K = zeros(2,2,N);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Knobs to turn
%--------------------------------------------------------------------------
% True initial state
x(:,1) = [ 0         ;
           3*pi/500 ];

% Initial aposteriori state estimate
x_aposteriori(:,1) = [ 1         ;
                       1*pi/500 ];

% Process noise covariance
Q = [   0.001    0 ;
          0      0 ];

% Measurement noise covariance
R = [ 0.1    0   ;
      0    0.01 ];

% Initial aposteriori error covariance estimate
P_aposteriori(:,:,1) = [ 1   0  ; 
                         0   1 ];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Perform extend Kalman filtering
%--------------------------------------------------------------------------
for i = 2:N
    % Update true state
    x(:,i) = [sin(x(2,i-1)*(i-1));x(2,i-1)] + sqrt(Q)*[randn;randn];

    % Update measurements
    z(:,i) = x(:,i) + sqrt(R)*[randn;randn];

    % Update apriori estimate
    x_apriori(:,i) = [sin(x_aposteriori(2,i-1)*(i-1));x_aposteriori(2,i-1)];

    % Update state Jacobian
    Ai = [0   (i-1)*cos(x_aposteriori(2,i-1)*(i-1)) ;
          0                    1                    ];

    % Assume knowledge of Q and R (Use system I.D. techniques in practice)
    Qi = Q;
    Ri = R;

    % Update aprioiri error covariance estimate
    P_apriori(:,:,i) = Ai*P_aposteriori(:,:,i-1)*Ai' + W*Qi*W';

    % Update Kalman gain
    K(:,:,i) = P_apriori(:,:,i)*H' / (H*P_apriori(:,:,i)*H'+V*Ri*V');

    % Update aposteriori state estimate
    x_aposteriori(:,i) = x_apriori(:,i) + K(:,:,i) * (z(:,i) - x_apriori(:,i));

    % Update aposteriori error covariance estimate
    P_aposteriori(:,:,i) = (eye(2) - K(:,:,i)*H) * P_apriori(:,:,i);
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------------
figure

subplot(2,1,1)
% Actual state position
b = plot(2:N,x(1,2:N),'b');
hold on
% State position estimates
r = plot(2:N,x_aposteriori(1,2:N),'r');
% State measurements
g = plot(2:N,z(1,2:N),'g+');
title('Extended Kalman Filtering of a Sine Wave - Position');
legend([b r g],'True Position','Position Estimates','Position Measurements');
xlabel('Time')
ylabel('Position')
grid on

subplot(2,1,2)
% Actual state frequency
b = plot(2:N,x(2,2:N),'b');
hold on
% State frequency estimates
r = plot(2:N,x_aposteriori(2,2:N),'r');
% Frequency measurements
g = plot(2:N,z(2,2:N),'g+');
title('Extended Kalman Filtering of a Sine Wave - Frequency');
legend([b r g],'True Frequency','Frequency Estimates','Frequency Measurements');
xlabel('Time')
ylabel('Frequency')
grid on
%--------------------------------------------------------------------------
