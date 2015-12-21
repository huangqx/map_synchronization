function [mid] = robust_median(x)

w = ones(1, length(x));
for iter = 1:1000
    mid = sum(w.*x)/sum(w);
    tp = x - mid;
    tp = (tp.*tp);
    ids = find(w > 0);
    [s,id] = max(tp(ids));
    if max(tp(ids)) < 1e-3
        break;
    end
        
    id = ids(id);
    w(id) = 0;
end
[x;w]
fprintf('mid = %f.\n', mid);
cvx_begin quiet
variable y(1)
minimize norm(x-y, 1)
cvx_end
fprintf('mid = %f.\n', y);