function [left_flats, ranks] = find_biggest_flattenings(dims, m)

if nargin < 2
    m = 1;
end
   
log_dims = log(dims);

n_dims = length(dims);
subset = [false(1, n_dims-1), true];

left_flats = true(0, n_dims);
k = 0;
ranks = [];

keep_running = true;
while keep_running
    
    if sum(subset)>1
        ld_left = log_dims(subset);
        lp_left = sum(ld_left);
        lp_left = lp_left + log1p(-sum(exp(ld_left-lp_left)));
        flat_rank = min(lp_left, sum(log_dims(~subset)));
        
        if k<m || flat_rank > ranks(m)
            if k<m
                k = k + 1;
            end
            i = search_sorted(ranks, flat_rank);
            ranks = [ranks(1:i-1), flat_rank, ranks(i:k-1)];
            left_flats = [left_flats(1:i-1,:); ...
                          subset; left_flats(i:k-1,:)];
    
        end
    end

    for i=n_dims:-1:1
        if subset(i) == 0
            subset(i) = true;
            subset(i+1:end) = false;
            break
        elseif i==1
            keep_running = false;
        end
    end
            
end

ranks = round(exp(ranks));

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

