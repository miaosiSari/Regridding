function y=w_vector(N, alpha)
if ~exist('alpha', 'var')
    alpha = 2.55;
end
y = zeros(N, 1);
for cnt = 1:N
    n = cnt-1;
    y(cnt) = w_single(n, N, alpha);
end
