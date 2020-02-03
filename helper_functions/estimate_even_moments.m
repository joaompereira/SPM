function [varargout] = estimate_even_moments(X)
% Estimate even moments from data
% Only even moments are needed to debias even moments

[L, N] = size(X);

maxnumel = 2e8;

n2 = nargout;

assert(n2==abs(round(n2-1))+1,'Are you even positive that n is an even positive integer?');

maxN = floor(maxnumel/L.^ceil(n2));

if maxN>N
    M_est = sum_tensor_power(X,n2);
else
    % avoid run-out-of-memory error
    K = floor(N/maxN);
    if N > maxN*K
        M_est = sum_tensor_power(X(:,maxN*K+1:N),n2);
    else
        M_est(1:n2) = {0};
    end
    for k=0:K-1
        M_est_ = sum_tensor_power(X(:,maxN*k+1:maxN*(k+1)),n2);
        for i = 1:n2
            M_est{i} = M_est{i} + M_est_{i};
        end
    end
end
   
for i = 1:n2
    M_est{i} = M_est{i}/N;    
end

varargout = M_est;

end

function M = sum_tensor_power(X, n2)

L = size(X,1);

M = cell(1,n2);

Xpow = X;
M{1} = X*X';

% Use kron to calculate tensor powers up to order n/2
for j=2:n2
    Xpow = khatri_rao_product(Xpow,X);
    % Speed up calculation using matrix multiplication
    M{j} = reshape(Xpow*Xpow',L*ones(1,2*j));
end

end


function AB = khatri_rao_product(A,B)
% Kronecker product of columns of 2 matrices
%  inspired on MATLAB kron
   n = size(A,2);

   A = reshape(A,1,[],n);
   B = reshape(B,[],1,n);
   AB = reshape(A.*B,[],n);

end    
