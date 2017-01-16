function [x_prior U_prior D_prior] = thornton(x_post,U_post,D_post,Uq,Dq,varargin)
%--------------------------------------------------------------------------
% Syntax: [x_prior U_prior D_prior]=thornton(x_post,U_post,D_post,Uq,Dq);
%         [x_prior U_prior D_prior]=thornton(x_post,U_post,D_post,Uq,Dq,A);
%
% Inputs:       x_post is the aposteriori state estimate
%               [U_post D_post] = myUD(P_post);
%               [Uq Dq] = myUD(Q);
%               A is the state transition matrix (can exclude when A is an
%               identity matrix)
%
% Outputs:      x_prior is the apriori state estimate
%               [U_prior D_prior] = myUD(P_prior);
%
% Description:  This function performs the Thornton square root Kalman
%               filter time update, which employs the modified weighted QR
%               decomposition. That is, it performs
%
%               x_prior = A * x_post;
%               P_prior = A * P_post * A' + Q;
%
%               but returns x_prior, U_prior, and D_prior, where
%               [U_prior D_prior] = myUD(P_prior);
%
% Author:       Brian Moore
%               brimoor@umich.edu
%
% Date:         June 28, 2012
%--------------------------------------------------------------------------

tol = 1e-15;

n = length(x_post);
sigma = 0; %#ok
dinv = 0; %#ok
i = 0; %#ok
j = 0; %#ok
k = 0; %#ok

a1 = zeros(1,n);
a2 = zeros(1,n);
v1 = zeros(1,n);
v2 = zeros(1,n);
D_prior = zeros(n);

if nargin == 6
    A = varargin{1};
    
    x_prior = A * x_post;

    % Form U_prior = A * U_post
    U_prior = A;
    for i = 1:n
        for j = n:-1:1
            sigma = U_prior(i,j);
            for k = 1:(j-1)
                sigma = sigma + U_prior(i,k) * U_post(k,j);
            end
            U_prior(i,j) = sigma;
        end
    end
else
    x_prior = x_post;
    U_prior = U_post;
end

for j = n:-1:1
    sigma = 0;
    for k = 1:n
        v1(k) = U_prior(j,k);
        a1(k) = D_post(k,k) * v1(k);
        sigma = sigma + v1(k) * a1(k);
    end
    for k = 1:n
        v2(k) = Uq(j,k);
        a2(k) = Dq(k,k) * v2(k);
        sigma = sigma + v2(k) * a2(k);
    end
    U_prior(j,j) = sigma;
    if sigma < tol
        error('New error covariance matrix is not positive definite');
    end
    dinv = 1 / sigma;
    for k = 1:(j-1)
        sigma = 0;
        for i = 1:n
            sigma = sigma + U_prior(k,i) * a1(i);
        end
        for i = 1:n
            sigma = sigma + Uq(k,i) * a2(i);
        end
        sigma = sigma * dinv;
        for i = 1:n
            U_prior(k,i) = U_prior(k,i) - sigma * v1(i);
        end
        for i = 1:n
            Uq(k,i) = Uq(k,i) - sigma * v2(i);
        end
        U_prior(j,k) = sigma;
    end
end

for j = 1:n
    D_prior(j,j) = U_prior(j,j);
    U_prior(j,j) = 1;
    for i = 1:(j-1)
        U_prior(i,j) = U_prior(j,i);
        U_prior(j,i) = 0;
    end
end
