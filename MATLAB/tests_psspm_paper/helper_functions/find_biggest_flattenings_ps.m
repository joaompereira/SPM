function [left_flats, ranks] = find_biggest_flattenings_ps(dims, nsyms, m)

if nargin < 2
    nsyms = ones(size(dims));
end
if nargin < 3
    m = 3;
end

assert(all(size(nsyms)==size(dims)))
   

n_dims = length(dims);
subset = [zeros(1, n_dims-1), 1];
left_flats = zeros(0, n_dims);
k = 0;
ranks = [];

while true
    left_rank = sv_dimension(dims, subset) - sum(dims(subset>0));
    right_rank = sv_dimension(dims, nsyms - subset);
    flat_rank = min(left_rank, right_rank);

    if k<m || flat_rank > ranks(m)
        if k<m
            k = k + 1;
        end
        i = search_sorted(ranks, flat_rank);
        ranks = [ranks(1:i-1), flat_rank, ranks(i:k-1)];
        left_flats = [left_flats(1:i-1,:); ...
                      subset; left_flats(i:k-1,:)];

    end
    
    for i=n_dims:-1:1
        if subset(i) < nsyms(i)
            subset(i) = subset(i) + 1;
            subset(i+1:end) = 0;
            break
        end
    end
    if all(subset==nsyms)
        break
    end
            
end

end


function [low] = search_sorted(v, val)
    
    n = length(v);

    low = 1;
    high = n+1;

    while high>low
        m = floor((low+high)/2);
        if val<v(m)
            low = m+1;
        else
            high = m;
        end
    end    


end



