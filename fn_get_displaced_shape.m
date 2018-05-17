function uu = fn_get_displaced_shape(xx, x, u, S, elements, k) 
%gets displaced shape at any point by back substitution into shape
%functions
const_per_el = 4;
ABCD = S * u;
uu = zeros(size(xx));
for ei = 1:size(elements, 1)
    ii = find(xx >= x(elements(ei, 1)) & xx <= x(elements(ei, 2)));
    jj = fn_el_to_global_index(ei, [1:const_per_el]);
    x0 = mean(x(elements(ei,:)));
    uu(ii) = fn_flex_shape_function(k(ei), xx(ii) - x0) * ABCD(jj);
end
end