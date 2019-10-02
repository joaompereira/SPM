function F = FOOBI1_tensor_decomposition(T, L, R, eigtol)
% Returns rank decomposition of order 4 symmetric tensors using FOOBI1

% Tolerance for choosing rank using eigenvalues
if ~exist('eigtol','var') || isempty(eigtol)
    eigtol = 1e-3;
end

if ~exist('L','var') || isempty(L)
    L = size(T,1);
end

% Flatten T and compute eigenvectors
% This is more stable than Cholesky decomposition if T is noisy
[V,D] = eig2(reshape(T,L^2,L^2));

if ~exist('R','var') || isempty(R)
    R = find(D>eigtol,1,'last');
end

H = V(:,1:R).*sqrt(D(1:R)');

% Get R(R+1)/2 unique indices of symmetric RxR matrix
ind = nchoosek(1:R+1,2)-[0,1];
indx = ind(:,1);
indy = ind(:,2);

R2 = nchoosek(R+1,2);

% To find the null space of P we calculate the eigenvectors of P'*P, and
% return the ones associated with smaller eigenvalues.
% To calculate P'*P we use the following formula
% (P'*P)_{ab,cd} = 4 * tr(H_a H_c) tr(H_b H_d) + 
%                  4 * tr(H_a H_d) tr(H_b H_c) - 
%                  8 * tr(H_a H_c H_b H_d)
% We also use that tr(H_a H_c) = 0 if a!=c, since the columns of H are
% orthogonal

HHdiag = D(indx).*D(indy)+(indx==indy).*D(indx).^2;
% Vectorized multiplication
% HH has all products of the form H_a H_c
HH = reshape(sum(reshape(H,L,L,1,R).*reshape(H,L,1,L,1,R)),L^2,R^2);
% The entries of HH2 are tr(H_a H_c H_b H_d) for different a,b,c,d
HH2 = HH'*HH;

% Relevant a,b,c,d indices
indices = indx + (indx'-1)*R + (indy'-1)*R^2 + (indy-1)*R^3;
% Ensuring PP is symmetric
indices = min(indices,indices');

PP = diag(HHdiag) - 2*HH2(indices);

% Using eigs here is faster than eig
% U is the null space of P
[U, ~] = eigs(PP+1e-5*eye(R2),R,'smallestabs');

W = zeros(R,R,R);

ind=0;

% Set W as the collection of matrices in the null space of P
for i=1:R
    ind = ind + 1;
    W(i,i,:) = 2*U(ind,:);
    for j=i+1:R
        ind = ind + 1;
        W(i,j,:) = U(ind,:);
        W(j,i,:) = U(ind,:);
    end
end

% Perform simultaneous diagonalization using Jacobi rotations
[Q, ~] = simultaneous_diagonalization(W);

A = H * Q;

F = zeros(L,R);

for i=1:R
    % Columns of f are rank 1 approximations of A
    MatA = reshape(A(:,i),L,L);
    [U, S, ~] = svd(MatA);
    
    F(:,i) = U(:,1)*sqrt(S(1,1));
end

end

function [Q, W] = simultaneous_diagonalization(W)
% Find simultaneous diagonaling orthogonal matrix Q with Jacobi iterations

% W is a set of K LxL matrices
[L, L2, K] = size(W);

tol = eps;

assert(L2 == L, 'Matrix not symmetric')

% Start with a good initialization for the orthogonal matrix
[Q, ~] = eig2(sum(W.*rand(1,1,K),3));

% Rotate matrices in W accordingly
for k=1:K
    W(:,:,k) = Q' * W(:,:,k) * Q;
end

s=1;

% Perform Jacobi iterations using formula in: JF Cardoso, A Souloumiac,
%  'Jacobi Angles For Simultaneous Diagonalization', 
%   Stop condition: the rotation angle is almost 0 (|sin x| = |s| < tol)
while abs(s) > tol
  
    off = triu(vecnorm(W,2,3),1);
    
    % calculate index of maximum off-diagonal element
    [offx, ind] = max(off);
    [~, j] = max(offx);
    i = ind(j);

    G = 0;
    
    % Calculate matrix for eigenvector calculation
    for k=1:K
        h = [W(i,i,k)-W(j,j,k); 2*W(i,j,k)];
        G = G + h*h';
    end

    [V,~] = eig(G,'vector');

    v =  V(:,2);
    % Ensure no abrupt changes in signal
    v = sign(v(1))*v/norm(v);

    % Calculate cosine (c) and sine (s) of Jacobi Rotation 
    c = sqrt(v(1)/2+.5);
    s = v(2)/sqrt(2*(v(1)+1));
    
    % R is the Jacobi rotation
    R = eye(L);
    R([i j],[i j]) = [c -s; s c];

    % Perform Jacobi rotation to matrices W
    for k=1:K
        W(:,:,k) = R' * W(:,:,k) * R;
    end

    Q = Q * R;

end

end
