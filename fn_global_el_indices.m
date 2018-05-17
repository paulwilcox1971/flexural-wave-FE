function ii = fn_global_el_indices(el_no)
const_per_el = 4;
ii = [(el_no - 1)* const_per_el + 1: el_no * const_per_el];
end

