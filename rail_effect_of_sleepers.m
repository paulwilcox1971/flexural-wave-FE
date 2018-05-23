clear;
close all;

%This script demonstrates using fn_calculate_emergent_waveno_for_model to 
%calculate the apparent dispersion curve for a rail with sleepers under 
%different loads, compared to the free rail under load case 

%Frequencies and tensions to model
freq = logspace(-2,2,100);
T = [-1,0,1] * 30e3*10;

%Rail properties
E = 210e9;
I = 30e-6;
m = 67.7;

col = 'rkb';%colours to plot different loads (needs to be same length as T)

%Sleeper details (is fastest and works fine with just one sleeper in model)
no_sleepers_in_model = 10;
sleeper_spacing = 0.65;
lateral_mass = 120;
lateral_stiffness = 50e6;
lateral_damping = 1e6;
rotational_mass = 0;
rotational_stiffness = 0;
rotational_damping = 0;

%--------------------------------------------------------------------------
%Define the model - first work out shortest wavelength to be encountered to
%determine max element size
[~, max_free_rail_waveno] = fn_waveguide_in_tension_dispersion(max(freq), E * I, min(T), m);
min_free_rail_wavelength = 2 * pi / max_free_rail_waveno;
min_nodes_per_wavelength = 1; %NB you don't need to change this as it doesn't alter accuracy - only reason why elements can't be too big is because of numerical instability due to exp(kx) terms
max_node_spacing = min_free_rail_wavelength / min_nodes_per_wavelength;
[nodes, elements, sleeper_nodes, forcing_node] = fn_create_rail_mesh(no_sleepers_in_model * sleeper_spacing, sleeper_spacing, [], max_node_spacing, 'ends at midspans');

%Run the modal calculation for all frequencies and tensions of interest
emergent_waveno = zeros(length(freq), length(T));
free_rail_waveno = zeros(length(freq), length(T));
for ti = 1:length(T)
    [~, free_rail_waveno(:, ti)] = fn_waveguide_in_tension_dispersion(freq, E * I, T(ti), m);
    for fi = 1:length(freq)
        clear('BC');
        for bi = 1:length(sleeper_nodes)
            BC(bi) = fn_BC_values_for_sleeper(lateral_mass, lateral_stiffness, lateral_damping, rotational_mass, rotational_stiffness, rotational_damping, freq(fi), sleeper_nodes(bi));
        end
        [K, S] = fn_build_flex_global_matrices(nodes, elements, ones(size(elements,1)) * E * I, ones(size(elements, 1), 1) * free_rail_waveno(fi, ti));
        emergent_waveno(fi, ti) = fn_calculate_emergent_waveno_for_model(K, S, BC, nodes, elements, ones(size(elements, 1), 1) * free_rail_waveno(fi, ti));
    end
end

%--------------------------------------------------------------------------
%Plot results
figure;
for ti = 1:length(T)
    h(ti) = loglog(freq, 2 * pi * freq(:) ./ real(emergent_waveno(:, ti)), [col(ti), '-']);
    hold on;
    loglog(freq, 2 * pi * freq(:) ./ real(free_rail_waveno(:, ti)), [col(ti), ':']);
end
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');
legend(h, num2str(T'))
return
figure;
for ti = 1:length(T)
    h(ti) = semilogx(freq, abs(imag(emergent_waveno(:, ti))), [col(ti), '-']);
    hold on;
end
xlabel('Frequency (Hz)');
ylabel('Attenuation (Np/m)');
legend(h, num2str(T'))