function [M] = generate_lowrank_tensor(X, n, lambda)

[L, R] = size(X);
n2 = floor(n/2);

assert(n==abs(round(n-1))+1,'Are you positive n is a positive integer?');

if nargin<3
  lambda = ones(1,R);
end

if n==1
  M = sum(X,2);
  return
end

Xpow = X;

% Use kron to calculate tensor powers up to order n/2
for j=2:n2
    Xpow = khatri_rao_product(Xpow,X);
end

% Speed up calculation using matrix multiplication
if mod(n,2)
    M = khatri_rao_product(Xpow,X)*(lambda.*Xpow)';
else
    M = Xpow*(lambda.*Xpow)';
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