function ii = fn_el_to_global_index(el_no, const_no)
%Converts element number and constant number (i.e. A=1,B=2,C=3,D=4) to
%indices into global S matrix.
%Either el_no or const_no can be vector, but not both.
const_per_el = 4;
ii = [(el_no - 1)* const_per_el + 1: el_no * const_per_el];
end

