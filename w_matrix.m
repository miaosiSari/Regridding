function y = w_matrix(N, alpha)

if ~exist('alpha', 'var')
    alpha = 2.55;
end

y1 = w_vector(N, alpha);
y2 = w_vector(N, alpha);
y = zeros(N);
for cnt1 = 1:N
    for cnt2 = 1:N
        y(cnt1, cnt2) = y1(cnt1) * y2(cnt2); 
    end
end