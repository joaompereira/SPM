function   [A,B,filler] = svd_ortho(T,r)
    % jennrich algorithm for partially symmetric tensor with symmetry in the first two indices 
    [~,n,k] = size(T);
    Tmatrix = reshape(T,[],k);
    Tmatrix1 = reshape(T,n,[]);
    [U,~,~] = svd(Tmatrix1);
    A = U(:,1:r);
    krp = reshape(A, [], 1, r) .* reshape(A, 1, [], r);
    krp = reshape(krp, [], r);
    Btranspose = pinv(krp)*Tmatrix;
    B = Btranspose';
    filler = [];
end