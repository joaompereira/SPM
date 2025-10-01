function [A,B,filler] = tensorlab_nls(T, r)
    n  = size(T,1);
    k = size(T,3);
    model=struct;
    model.variables.u=randn(n,r);
    model.variables.v=randn(k,r);
    model.factors.U='u';
    model.factors.V='v';
    model.factorizations.myfac.data = T;
    model.factorizations.myfac.cpd={'U','U','V'};
    sol = ccpd_nls(model);
    A = sol{1};
    B = sol{2}.* vecnorm(A).^2;
    A = A./vecnorm(A);
    filler = [];
end