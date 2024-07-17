function [X] = tensorlab_ccpd_nls(T, L, n, R, U0, maxiter)
% Rank decomposition using TensorLAB
% TensorLab package needed to run this algorithm
% You can get it at 'www.tensorlab.net'
% We use ccpd_nls since, empirically, it converged faster

% NOTE: Tensorlab has their own implementation of the Matlab routines
% kron and kmeans, which are used in SPM or GPCA. It is recommended to
% rename these functions in the Tensorlab folder so there is no conflict 
% between the Matlab and the Tensorlab implementations.

if nargin < 2 || isempty(L)
    L = size(T,1);
end

if nargin < 3 || isempty(n)
    n = round(log(numel(T))/log(L));
end

if nargin < 4|| isempty(R)
    error('TensorLab cannot estimate rank of tensor')
end

if nargin < 5 || isempty(U0)
    U0 = randn(L,R);
end

if nargin < 6 || isempty(maxiter)
    maxiter = 500;
end


% Definitions for TensorLab ccpd_nls
model = struct;
model.variables = U0; % Rank R Tensor in L variables
model.factors = 1;
model.factorizations.unconst.data = T;
model.factorizations.unconst.cpd  = ones(1,n); % Tensor order

% Some definitions to ensure same conditions as other algorithms.
% Since this is a second order method, once it finds the global minimum
% the logarithm of the error decreases quadratically
Ures = ccpd_nls(model,'MaxIter', maxiter,...
                      'TolFun',1e-24,...
                      'TolX', 1e-12);

X = Ures{1};

