function [err] = rderror_normalized(A_est, A_true, varargin)
% Returns the error of the X_est in terms of X_true
    
    if nargin>=5
        lambda_est = varargin{1};
        lambda_true = varargin{2};
        n = varargin{3};
    else
        lambda_est = ones(1, size(A_est, 2));
        lambda_true = ones(1, size(A_true, 2));
        n = varargin{1};
    end
    
    T = generate_lowrank_tensor(A_true, lambda_true, n);
    T_est = generate_lowrank_tensor(A_est, lambda_est, n);
 
    err = norm(T(:) - T_est(:)) / norm(T(:));
    
end