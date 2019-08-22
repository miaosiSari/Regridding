function y = w_single(n, N, alpha)

if ~exist('alpha', 'var')
    alpha = 2.55;
end
numerator = pi*alpha*sqrt(1-(2*n/(N-1)-1)^2);
denominator = pi*alpha;
y = modified_Bessel(numerator)/modified_Bessel(denominator);