function [sigma] = estimate_sigma_4thmoment(M2, M4)

L = size(M2,1);
d = L*(L+1)/2;

[V, d4] = eig2(reshape(M4,L^2,L^2));

vI = reshape(eye(L),1,[]);

IM2 = symmetrize_tensor(M2(:)*vI,L,4);

I4 = symmetrize_tensor(vI'*vI,L,4);

a = V(:,d)'*(I4*V(:,d));

b = V(:,d)'*(IM2*V(:,d));

c = d4(d);

sigma = sqrt((b-sqrt(b^2-a*c/3))/a);