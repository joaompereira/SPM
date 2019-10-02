function [M] = generate_lowrank_tensor(X, n)

[L, R] = size(X);
n2 = floor(n/2);

assert(n==abs(round(n-1))+1,'Are you positive n is a positive integer?');

if n==1
  M = sum(X,2);
  return
end

Xpow = X;

% Use kron to calculate tensor powers up to order n/2
for j=2:n2
    Xpow = columnkron(Xpow,X);
end

% Speed up calculation using matrix multiplication
if mod(n,2)
    M = columnkron(Xpow,X)*Xpow';
else
    M = Xpow*Xpow';
end

M = reshape(M,L*ones(1,n));

end


function AB = columnkron(A,B)
% Kronecker product of columns of 2 matrices
%  inspired on MATLAB kron
   [a,n] = size(A);
   b = size(B,1);
   if size(B,2)~=n
       error('Different number of columns between A and B!');
   end

   A = reshape(A,[1 a n]);
   B = reshape(B,[b 1 n]);
   AB = reshape(A.*B,[a*b n]);

end