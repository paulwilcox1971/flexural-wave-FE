function [xx, uu] = fn_get_displaced_shape(nodes, elements, u, S, k, pts_per_element) 
%gets displaced shape at any point by back substitution into shape
%functions
min_element = min(abs(nodes(2:end) - nodes(1:end-1)));
xx = linspace(min(nodes), max(nodes), (max(nodes) - min(nodes)) / (min_element / pts_per_element));



const_per_el = 4;
ABCD = S * u;
uu = zeros(size(xx));
for ei = 1:size(elements, 1)
    ii = find(xx >= nodes(elements(ei, 1)) & xx <= nodes(elements(ei, 2)));
    jj = fn_el_to_global_index(ei, [1:const_per_el]);
    x0 = mean(nodes(elements(ei,:)));
    uu(ii) = fn_flex_shape_function(k(ei), xx(ii) - x0) * ABCD(jj);
end
end