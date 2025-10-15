function [lambda,factors_norm] = tensorlab_gevd_111(T, r)
    sol = cpd_gevd(T,r);
    A = sol{1};
    B = sol{2};
    C = sol{3};
    A_ = A./vecnorm(A);
    B_ = B./vecnorm(B);
    C_ = C./vecnorm(C);
    lambda = vecnorm(A).*vecnorm(B).*vecnorm(C);
    factors_norm = {A_,B_,C_};
end