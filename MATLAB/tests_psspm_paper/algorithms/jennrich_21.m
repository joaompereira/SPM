function   [A,B,filler] = jennrich_21(T,r)
    % jennrich algorithm for partially symmetric tensor with symmetry in the first two indices 
    [~,~,k] = size(T);
    z = randn(k,1);
    zp = randn(k,1);
    Mz   = sum(bsxfun(@times, T, reshape(z,1,1,[])), 3);
    Mzp  = sum(bsxfun(@times, T, reshape(zp,1,1,[])), 3);
    [A,D] = eig(Mz*pinv(Mzp));
    [~, idx] = sort(abs(diag(real(D))), 'descend');
    A = A(:,idx(1:r));
    A = real(A)./vecnorm(real(A));
    krp = reshape(A, [], 1, r) .* reshape(A, 1, [], r);
    krp = reshape(krp, [], r);
    Tmatrix = reshape(T,[],k);
    Btranspose = pinv(krp)*Tmatrix;
    B = Btranspose';
    filler = [];
end

