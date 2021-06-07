function [err] = rderror_general(A_est, A_true, W_est, W_true, n)

 T = generate_lowrank_tensor_general(A_true, W_true, n);
 T_est = generate_lowrank_tensor_general(A_est, W_est, n);
 
 err = norm(T(:) - T_est(:));
 
end