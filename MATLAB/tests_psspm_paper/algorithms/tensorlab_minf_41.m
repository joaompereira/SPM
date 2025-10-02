function [lambda,factors_norm] = tensorlab_minf_41(T, r)
    n  = size(T,1);
    k = size(T,5);
    model=struct;
    model.variables.u=randn(n,r);
    model.variables.v=randn(k,r);
    model.factors.U='u';
    model.factors.V='v';
    model.factorizations.myfac.data = T;
    model.factorizations.myfac.cpd={'U','U','U','U','V'};
    sol = ccpd_minf(model);
    A = sol{1};
    B = sol{2};
    A_ = A./vecnorm(A);
    B_ = B./vecnorm(B);
    lambda = vecnorm(A).^4.*vecnorm(B);
    factors_norm = {A_,B_};
end