function [A,B,filler] = tensorlab_als(T, r)
    U = randn(size(T,1), r);
    V = randn(size(T,3), r);
    sol = cpd_als(T,{U,U,V});
    A = sol{1};
    A = A./vecnorm(A);
    filler = [];
    krp = reshape(A, [], 1, r) .* reshape(A, 1, [], r);
    krp = reshape(krp, [], r);
    k = size(T,3);
    Tmatrix = reshape(T,[],k);
    Btranspose = pinv(krp)*Tmatrix;
    B = Btranspose';
end