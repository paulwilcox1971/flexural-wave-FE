function v = fn_damper_values_for_reflection_coefficient(RC, EI, k)
v = [k ^ 3 * EI / (1 - RC), k * EI * (1 - RC)];
end