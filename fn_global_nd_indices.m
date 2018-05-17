function ii = fn_global_nd_indices(nd_no)
dof_per_nd = 2;
ii = [(nd_no - 1)* dof_per_nd + 1: nd_no * dof_per_nd];
end
