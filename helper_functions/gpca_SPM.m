function [groups, S] = gpca_SPM(X, k, n)

[L, K] = size(X);
d  = 4;
d2 = d/2;

[M{1:d2}] = estimate_even_moments(X);

if d==4
    sigma = estimate_sigma_4thmoment(M{1}, M{2});
    [dM{1:d2}] = debias_even_moments(L, sigma, M{:});
else
    dM = M;
end

tries = 5;

error = Inf;

for i = 1:tries

    S = subspace_power_method_general(dM{d2}, L, d, n*nchoosek(k+d2-1,d2), k);

    n = length(S);
    dist2S = zeros(n, K);
    for i=1:n
      dist2S(i,:) = vecnorm(X - S{i}*S{i}'*X);
    end

    [mdist, groups_] = min(dist2S);

    if norm(mdist)<error
       groups = groups_;
       error = norm(mdist);
    end

end


end