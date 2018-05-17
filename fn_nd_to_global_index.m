function ii = fn_nd_to_global_index(nd, dof)
%Converts node and DOF (1 or 2) to index into global vectors and matrices.
%Either node or DOF can be a vector, but not both
dof_per_nd = 2;
ii = (nd - 1) * dof_per_nd + dof;
end