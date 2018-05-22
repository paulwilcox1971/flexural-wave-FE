function BC = fn_BC_values_for_lateral_mass_loading(mass, freq, node)
BC.node = node;
BC.value = -mass * (2 * pi * freq) ^ 2;
BC.type = 'lateral stiffness';
end