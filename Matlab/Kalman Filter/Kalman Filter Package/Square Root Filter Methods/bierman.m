function [x_post U_post D_post] = bierman(z,R,H,x_prior,U_prior,D_prior)
%--------------------------------------------------------------------------
% Syntax:  [x_post U_post D_post] = bierman(z,R,H,x_prior,U_prior,D_prior);
%
% Inputs:       z is a scalar measurement
%               R is the variance of z
%               H is a row vector of length length(x_prior)
%               x_prior is the apriori state estimate
%               [U_prior D_prior] = myUD(P_prior);
%               
% Outputs:      x_post is the aposteriori state estimate
%               [U_post D_post] = myUD(P_post);
%
% Description:  This function performs the Bierman square root Kalman
%               filter scalar measurement update. That is, it performs
%
%               K = P_prior * H' / (H * P_prior * H' + R);
%               x_post = x_prior + K * (z - H * x_prior);
%               P_post = (I - K*H) * P_prior;
%
%               but returns x_post, U_post, and D_post, where
%               [U_post D_post] = myUD(P_post);
%
% Author:       Brian Moore
%               brimoor@umich.edu
%
% Date:         June 28, 2012
%--------------------------------------------------------------------------

x_post = x_prior;
U_post = U_prior; 
D_post = D_prior;
a = U_post' * H';
b = D_post * a;
dz = z - H * x_prior; 
alpha = R; 
gamma = 1 / alpha; 
  for j = 1:length(x_prior)
  beta   = alpha; 
  alpha  = alpha + a(j) * b(j); 
  lambda = -a(j) * gamma; 
  gamma  = 1 / alpha; 
  D_post(j,j) = beta * gamma * D_post(j,j); 
    for i = 1:j-1 
    beta = U_post(i,j); 
    U_post(i,j) = beta + b(i) * lambda; 
    b(i) = b(i) + b(j) * beta; 
    end
  end
x_post = x_post + gamma * dz * b;
