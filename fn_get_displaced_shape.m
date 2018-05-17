function [y, y_at_nodes] = fn_get_displaced_shape(x, nodes, elements, ABCD, k) 

y = zeros(size(x));
y_at_nodes = zeros(size(nodes));
for ei = 1:size(elements, 1)
    ii = find(x >= nodes(elements(ei, 1)) & x <= nodes(elements(ei, 2)));
    jj = fn_global_el_indices(ei);
    y(ii) = fn_flex_shape_function(k(ei), x(ii)) * ABCD(jj);
    y_at_nodes(elements(ei,:)) = fn_flex_shape_function(k(ei), nodes(elements(ei,:))) * ABCD(jj);
end
