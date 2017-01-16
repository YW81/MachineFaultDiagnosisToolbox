%--------------------------------------------------------------------------
% Name:            KF_Demo2.m
%
% Description:     Quick examples of each of the four implementations of
%                  Kalman filters in this package:
%
%                  1. Standard Kalman Filter
%                  2. Square Root Kalman Filter
%                  3. Standard AR Kalman Filter
%                  4. Square Root AR Kalman Filter
%
% Author:          Brian Moore
%                  brimoor@umich.edu
%
% Date:            June 1, 2012
%--------------------------------------------------------------------------

%% 1. Standard Kalman Filter (Compare to #1 - should be identical!)

m = 3; % dimension
n = 300; % number of measurement samples
sigma = 1; % noise scaling
MAlen = 30; % moving average filter length
N = 15; % lookback window length for covariance estimation

% Fix the random number seed so we can compare performance of standard and
% square root implementations
rng(1);

% Generate true state
xmin = -4;
xmax = 3;
x = linspace(xmin,xmax,n);
y = randn(m,3) * [x ; x.^2 ; x.^3];

% Generate noisy measurements (with correlated Gaussian noise)
S = randn(m);
z = y + sigma * chol(S * S') * randn(m,n);

% Perform Standard Kalman filtering
[x_kf KF] = StandardKalmanFilter(z,MAlen,N);

% Plot results
figure
for i = 1:m
    subplot(m,1,i)
    p1 = plot(x((MAlen+N-1):end),z(i,(MAlen+N-1):end),'r+');
    hold on;
    grid on;
    p2 = plot(x((MAlen+N-1):end),y(i,(MAlen+N-1):end),'g');
    p3 = plot(x((MAlen+N-1):end),x_kf(i,(MAlen+N-1):end),'b');
    if (i == 1)
        title(['Standard Kalman Filter Example (' num2str(m) '-D Data)']);
    end
    if (i == ceil(m/2))
        legend([p1 p2 p3],'Noisy Measurements','True State','KF State Estimate');
    end
end

%% 2. Square Root Kalman Filter (Compare to #1 - should be identical!)

m = 3; % dimension
n = 300; % number of measurement samples
sigma = 1; % noise scaling
MAlen = 30; % moving average filter length
N = 15; % lookback window length for covariance estimation

% Fix the random number seed so we can compare performance of standard and
% square root implementations
rng(1);

% Generate true state
xmin = -4;
xmax = 3;
x = linspace(xmin,xmax,n);
y = randn(m,3) * [x ; x.^2 ; x.^3];

% Generate noisy measurements (with correlated Gaussian noise)
S = randn(m);
z = y + sigma * chol(S * S') * randn(m,n);

% Perform Square Root Kalman filtering
[x_kf KF] = SquareRootKalmanFilter(z,MAlen,N);

% Plot results
figure
for i = 1:m
    subplot(m,1,i)
    p1 = plot(x((MAlen+N-1):end),z(i,(MAlen+N-1):end),'r+');
    hold on;
    grid on;
    p2 = plot(x((MAlen+N-1):end),y(i,(MAlen+N-1):end),'g');
    p3 = plot(x((MAlen+N-1):end),x_kf(i,(MAlen+N-1):end),'b');
    if (i == 1)
        title(['Square Root Kalman Filter Example (' num2str(m) '-D Data)']);
    end
    if (i == ceil(m/2))
        legend([p1 p2 p3],'Noisy Measurements','True State','KF State Estimate');
    end
end

%% 3. Standard AR Kalman Filter (Compare to #4 - should be identical!)

m = 1; % dimension
n = 300; % number of measurement samples
sigma = 2; % noise scaling
MAlen = 30; % moving average filter length
N = 15; % lookback window length for covariance estimation
M = 2; % auto-regressive model order

% Fix the random number seed so we can compare performance of standard and
% square root implementations
rng(3);

% Generate true state
xmin = -4;
xmax = 3;
x = linspace(xmin,xmax,n);
y = randn(m,3) * [x ; x.^2 ; x.^3];
r = rand;
y(2:end) = r * y(1:(end-1)) + (1 - r) * y(2:end); % Add AR 

% Generate noisy measurements (with correlated Gaussian noise)
S = randn(m);
z = y + sigma * chol(S * S') * randn(m,n);

% Perform Standard AR Kalman filtering
[x_kf KF] = StandardARKalmanFilter(z,M,MAlen,N);

% Plot results
figure
for i = 1:m
    subplot(m,1,i)
    p1 = plot(x((MAlen+N-1):end),z(i,(MAlen+N-1):end),'r+');
    hold on;
    grid on;
    p2 = plot(x((MAlen+N-1):end),y(i,(MAlen+N-1):end),'g');
    p3 = plot(x((MAlen+N-1):end),x_kf(i,(MAlen+N-1):end),'b');
    if (i == 1)
        title(['Standard AR Kalman Filter Example (' num2str(m) '-D Data)']);
    end
    if (i == ceil(m/2))
        legend([p1 p2 p3],'Noisy Measurements','True State','KF State Estimate');
    end
end

%% 4. Square Root AR Kalman Filter (Compare to #3 - should be identical!)

m = 1; % dimension
n = 300; % number of measurement samples
sigma = 2; % noise scaling
MAlen = 30; % moving average filter length
N = 15; % lookback window length for covariance estimation
M = 2; % auto-regressive model order

% Fix the random number seed so we can compare performance of standard and
% square root implementations
rng(3);

% Generate true state
xmin = -4;
xmax = 3;
x = linspace(xmin,xmax,n);
y = randn(m,3) * [x ; x.^2 ; x.^3];
r = rand;
y(2:end) = r * y(1:(end-1)) + (1 - r) * y(2:end); % Add AR 

% Generate noisy measurements (with correlated Gaussian noise)
S = randn(m);
z = y + sigma * chol(S * S') * randn(m,n);

% Perform Square Root AR Kalman filtering
[x_kf KF] = SquareRootARKalmanFilter(z,M,MAlen,N);

% Plot results
figure
for i = 1:m
    subplot(m,1,i)
    p1 = plot(x((MAlen+N-1):end),z(i,(MAlen+N-1):end),'r+');
    hold on;
    grid on;
    p2 = plot(x((MAlen+N-1):end),y(i,(MAlen+N-1):end),'g');
    p3 = plot(x((MAlen+N-1):end),x_kf(i,(MAlen+N-1):end),'b');
    if (i == 1)
        title(['Square Root AR Kalman Filter Example (' num2str(m) '-D Data)']);
    end
    if (i == ceil(m/2))
        legend([p1 p2 p3],'Noisy Measurements','True State','KF State Estimate');
    end
end
