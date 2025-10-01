function sv_dim = sv_dimension(dims, nsyms)
% Calculate Segre-Veronese variety dimension
    sv_dim = sum(gammaln(dims + nsyms) - gammaln(dims) - gammaln(nsyms + 1));
    sv_dim = round(exp(sv_dim));
    
end