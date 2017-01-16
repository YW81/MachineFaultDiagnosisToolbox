function [x_aposteriori P_aposteriori] = KalmanFilterIteration(z,Q,R,x_aposteriori_last,P_aposteriori_last)
%
% Performs one iteration of the standard Kalman filter update equations
% Note: this function is used ONLY by KalmanFilterDemo.m
%

% Update apriori estimate
x_apriori = x_aposteriori_last;

% Update apriori error covariance estimate
P_apriori = P_aposteriori_last + Q;

% Update Kalman gain
% Note: backslash denotes right matrix inversion in MATLAB
K = P_apriori / (P_apriori + R);

% Update aposteriori state estimate
x_aposteriori = x_apriori + K * (z - x_apriori);

% Update aposteriori error covariance estimate
P_aposteriori = (eye(length(x_aposteriori)) - K) * P_apriori;
