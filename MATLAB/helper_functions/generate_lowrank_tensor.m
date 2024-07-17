function [M] = generate_lowrank_tensor(X, varargin)

[L, R] = size(X);

if nargin<3
  lambda = ones(1, R);
  n = varargin{1};
else
  lambda = varargin{1};
  n = varargin{2};
end

n2 = floor(n/2);

assert(n==abs(round(n-1))+1,'Are you positive n is a positive integer?');

if n==1
  M = sum(X,2);
  return
end

Xpow = khatri_rao_power(X, n2);

% Speed up calculation using matrix multiplication
if mod(n,2)
    M = khatri_rao_product(Xpow,X)*(lambda.*Xpow).';
else
    M = Xpow*(lambda.*Xpow).';
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

function An = khatri_rao_power(A, n)
    if n==1
        An = A;
    elseif n==2 
        An =  khatri_rao_product(A, A);
    else
        n2 = floor(n/2);
        An2 = khatri_rao_power(A, n2);
        if mod(n,2)
            An2_ = khatri_rao_product(An2, A);
        else
            An2_ = An2;
        end
        An = khatri_rao_product(An2, An2_);
    end

end