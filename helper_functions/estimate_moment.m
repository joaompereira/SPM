function [M_est] = estimate_moment(X, n)

[L, N] = size(X);

maxnumel = 2e8;

assert(n==abs(round(n-1))+1,'Are you positive n is a positive integer?');

maxN = floor(maxnumel/L.^ceil(n/2));

if n==1
    M_est = mean(X,2);
elseif maxN>N
    M_est = sum_tensor_power(X,n)/N;
else
    % avoid run-out-of-memory error
    K = floor(N/maxN);
    if N > maxN*K
        M_est = sum_tensor_power(X(:,maxN*K+1:N),n);
    else
        M_est = 0;
    end
    for k=0:K-1
        M_est = M_est + sum_tensor_power(X(:,maxN*k+1:maxN*(k+1)),n);
    end
    M_est = M_est/N;    
end

end

function M = sum_tensor_power(X, n)

L = size(X,1);
n2 = floor(n/2);

Xpow = X;

% Use kron to calculate tensor powers up to order n/2
for j=2:n2
    Xpow = khatri_rao_product(Xpow,X);
end

% Speed up calculation using matrix multiplication
if mod(n,2)
    M = khatri_rao_product(Xpow,X)*Xpow';
else
    M = Xpow*Xpow';
end

M = reshape(M,L*ones(1,n));

end


function AB = khatri_rao_product(A,B)
% Kronecker product of columns of 2 matrices
%  inspired on MATLAB kron
   n = size(A,2);

   A = reshape(A,1,[],n);
   B = reshape(B,[],1,n);
   AB = reshape(A.*B,[],n);

end    
