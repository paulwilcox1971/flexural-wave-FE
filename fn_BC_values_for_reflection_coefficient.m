function BC = fn_BC_values_for_reflection_coefficient(RC, EI, k, node)
BC.type = 'combined stiffness';
BC.node = node;
BC.value = 1i * [k ^ 3 * EI / (1 - RC), k * EI * (1 - RC)];
end