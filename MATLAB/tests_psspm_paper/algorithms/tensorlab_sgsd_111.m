function [lambda,factors_norm,err] = tensorlab_sgsd_111(T, r)
    U = randn(size(T,1), r);
    V = randn(size(T,2), r);
    W = randn(size(T,3),r);
    sol = cpd3_sgsd(T,{U,V,W});
    A = sol{1};
    B = sol{2};
    C = sol{3};
    A_ = A./vecnorm(A);
    B_ = B./vecnorm(B);
    C_ = C./vecnorm(C);
    factors = {A,B,C};
    T_est = generate_lowrank_tensor(ones(1,double(r)),factors{:}, [1,1,1]);
    err = norm(T-T_est,'fro');
    lambda = vecnorm(A).*vecnorm(B).*vecnorm(C);
    factors_norm = {A_,B_,C_};
end