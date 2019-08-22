function y = modified_Bessel(x, m)
%https://en.wikipedia.org/wiki/Bessel_function#Modified_Bessel_functions:_I%CE%B1,_K%CE%B1
if ~exist('m', 'var')
    m = 7;
end
y = 0.0;
for cnt = 0:m
    f = factorial(cnt);
    y = y + ((x / 2)^(2 * cnt))/(f^2);
end

