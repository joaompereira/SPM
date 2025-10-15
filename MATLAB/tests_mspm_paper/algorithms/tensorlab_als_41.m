function [lambda,factors_norm] = tensorlab_als_41(T, r)
    U = randn(size(T,1), r);
    V = randn(size(T,5), r);
    sol = cpd_als(T,{U,U,U,U,V});
    A = sol{1};
    B = sol{5};
    A_ = A./vecnorm(A);
    B_ = B./vecnorm(B);
    lambda = vecnorm(A).^4.*vecnorm(B);
    factors_norm = {A_,B_};
end