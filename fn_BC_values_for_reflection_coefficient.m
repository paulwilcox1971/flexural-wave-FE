function BC = fn_BC_values_for_reflection_coefficient(RC, EI, waveno_or_mode_shape, node)
BC.type = 'combined stiffness';
BC.node = node;
if length(waveno_or_mode_shape) == 1
    BC.value = 1i * [waveno_or_mode_shape ^ 3 * EI / (1 - RC), waveno_or_mode_shape * EI * (1 - RC)];
elseif length(waveno_or_mode_shape) == 4
    BC.value = [waveno_or_mode_shape(3) / waveno_or_mode_shape(1), waveno_or_mode_shape(4) / waveno_or_mode_shape(2)];
else
    error('Either wavenumber or 4-element mode shape needed for this BC');
end
end