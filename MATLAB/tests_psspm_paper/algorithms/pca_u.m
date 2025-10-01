function   [A,B,filler] = pca_u(T,r)
    n = size(T,1);
    k = size(T,3);
    mean_cov = zeros(n,n);
    for i = 1:k
        mean_cov = mean_cov + 1/k*T(:,:,i);
    end
    [A,~] = eig(mean_cov);
    A = A(:,1:r);
    krp = reshape(A, [], 1, r) .* reshape(A, 1, [], r);
    size(krp)
    Tmatrix = reshape(T,[],k);
    krp = reshape(krp, [], r);
    Btranspose = pinv(krp)*Tmatrix;
    B = Btranspose';
    filler = [];
end

    
    