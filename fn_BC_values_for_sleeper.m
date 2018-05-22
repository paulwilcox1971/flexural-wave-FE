function BC = fn_BC_values_for_sleeper(lateral_mass, lateral_stiffness, lateral_damping, rotational_mass, rotational_stiffness, rotational_damping, freq, node)
w = 2 * pi * freq;
BC.node = node;
BC.value = [lateral_stiffness + 1i * lateral_damping * w - lateral_mass * w ^ 2, ...
    rotational_stiffness + 1i * rotational_damping * w - rotational_mass * w ^ 2];
BC.type = 'combined stiffness';
end