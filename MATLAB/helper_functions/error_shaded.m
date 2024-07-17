function [hp, hf] = error_shaded(x, y, varargin)

x = x(:);

[b, a] = size(y);

assert(b == length(x), 'y wrong format');

if nargin<4 || ~isnumeric(varargin{2})
    assert(size(varargin{1},1) == b, 'err wrong format');
    assert(a == size(varargin{1},2), 'y and err diff num rows');
    err = varargin{1};
    err_up = y + err;
    err_down = y - err;
    varargin = varargin(2:end);
else
    err_up = varargin{1};
    err_down = varargin{2};
    for arg=varargin(1:2)
        assert(size(arg{1},1) == b, 'err wrong format');
        assert(a == size(arg{1},2), 'y and err diff num rows');
    end
    varargin = varargin(3:end);
end

hp = plot(x, y, varargin{:});
hold on
for i = a:-1:1
    err_upi = err_up(:, i);
    err_downi = err_down(:, i);
    nonanind = ~isnan(err_downi) & ~isnan(err_upi);
    xi = x(nonanind);
    err_upi = err_upi(nonanind);
    err_downi = err_downi(nonanind);

    hf(i) = patch([xi; xi(end:-1:1)]', [err_upi;...
            err_downi(end:-1:1)] , hp(i).Color);
    hf(i).EdgeColor = 'none';
    hf(i).FaceAlpha = 0.3;
end
hold off

set(gca, 'Children', flipud(get(gca, 'Children')) )

end