function [X] = tensorlab_ccpd_nls(T, L, n, R)
% Rank decomposition using TensorLAB
% TensorLab package needed to run this algorithm
% You can get it at 'www.tensorlab.net'
% We use ccpd_nls since, empirically, it converged faster

if ~exist('L','var') || isempty(L)
    L = size(T,1);
end

if ~exist('n','var') || isempty(n)
    n = round(log(numel(T))/log(L));
end

if ~exist('R','var') || isempty(R)
    error('TensorLab cannot estimate rank of tensor')
end

% Definitions for TensorLab ccpd_nls
model = struct;
% Rank R Tensor in L variables
model.variables = randn(L,R);
model.factors = 1;
model.factorizations.unconst.data = T;
% Tensor order
model.factorizations.unconst.cpd  = ones(1,n);

% Some definitions to ensure same conditions as other algorithms.
% Since this is a second order method, once it finds the global minimum
% the logarithm of the error decreases quadratically
Ures = ccpd_nls(model,'MaxIter',400,...
                      'TolFun',1e-24,...
                      'TolX', 1e-12);

X = Ures{1};

