function n = fn_flex_shape_function(k, x)
%returns basic shape function, n(x) terms for flexural element
%actual shape is y(x) = n(x) .* [A B C D] 
n = [exp(1i * k * x(:)), exp(-1i * k * x(:)), exp(k * x(:)), exp(-k * x(:))];
end