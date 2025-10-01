function [lambda,factors_norm,err] = tensorlab_als_41(T, r)
    U = randn(size(T,1), r);
    V = randn(size(T,5), r);
    sol = cpd_als(T,{U,U,U,U,V});
    A = sol{1};
    B = sol{5};
    A_ = A./vecnorm(A);
    B_ = B./vecnorm(B);
    factors = {A,B};
    T_est = generate_lowrank_tensor(ones(1,r),factors{:}, [4,1]);
    err = norm(T-T_est,'fro');

    lambda = vecnorm(A).^4.*vecnorm(B);
    factors_norm = {A_,B_};
end