function [left_flats, max_rank] = find_best_flatpair(dims, n_flats)

if nargin < 2
    n_flats = 3;
end

[left_flats, ranks] = find_biggest_flattenings(dims, n_flats);

for j=2:n_flats
    for i = 1:j-1
        if ~all(xor(left_flats(i, :),left_flats(j, :)))
            left_flats = left_flats([i, j], :);
            max_rank = ranks(j);
            return
        end
    end
end

error('Did not find a good flatpair');
   
end

