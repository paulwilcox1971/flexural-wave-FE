function [K, S] = fn_build_flex_global_matrices(nodes, elements, EI, k)
no_els = size(elements, 1);
dof_per_nd = 2;
nds_per_el = 2;
dof_per_el = nds_per_el * dof_per_nd;
const_per_el = 4;
no_nds = size(nodes, 1);

%get element matrices
Ke = zeros(dof_per_el, dof_per_el, no_els);
Se = zeros(dof_per_el, dof_per_el, no_els);
for ei = 1:no_els
    [Ke(:, :, ei), Se(:, :, ei)] = fn_flex_element_matrices([nodes(elements(ei,:))], EI(ei), k(ei));
end

%build global matrices
K = zeros(no_nds * dof_per_nd);
S = zeros(const_per_el * no_els, no_nds * dof_per_nd);
for ei = 1:no_els
    jj = [fn_nd_to_global_index(elements(ei, 1), [1, 2]), fn_nd_to_global_index(elements(ei, 2), [1, 2])];
    K(jj, jj) = K(jj, jj) + Ke(:, :, ei);
    ii = fn_el_to_global_index(ei, [1:const_per_el]);
    S(ii, jj) = Se(:, :, ei);
end

end