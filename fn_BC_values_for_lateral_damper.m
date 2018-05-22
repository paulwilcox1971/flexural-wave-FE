function BC = fn_BC_values_for_lateral_damper(damping, freq, node)
BC.node = node;
BC.value = -1i * damping * (2 * pi * freq);
BC.type = 'lateral stiffness';
end