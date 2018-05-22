clear;
close all;

%this script uses fn_calculate_emergent_waveno_for_model to calculate what
%the apparent dispersion curve is for a rail with point mass loading, M, representing
%sleepers. There is no damping or stiffness representing the sleeper
%supports so the sleepers just provide extra inertial loading = -w^2*M*u at
%each sleeper position; if sleeper support stiffness was included, there
%would be a +k*u term counteracting this. It would be easy to include a damping term, C, as well with
%another -i*w*C*u

%some rail-like parameters - calculate wavenumber in free rail under
%various tensions
freq = logspace(-1,2,20);
E = 210e9;
I = 6e-6;
m = 67.7;

%tensions too consider
T = [-1,0,1] * 30e3*10;
col = 'rkb';%colours to plot

%length to consider, sleepers etc
no_sleepers_in_model = 20;
sleeper_spacing = 0.6;
lateral_mass = 10;
lateral_stiffness = 0;
lateral_damping = 0.0001;
rotational_mass = 0;
rotational_stiffness = 0;
rotational_damping = 0;


max_node_spacing = 0.25;

%--------------------------------------------------------------------------
%make the track, one node per sleeper, with extra half sleeper pitch nodes
%at either end
[nodes, elements, sleeper_nodes, forcing_node] = fn_create_rail_mesh(no_sleepers_in_model * sleeper_spacing, sleeper_spacing, [], max_node_spacing, 'ends at midspans');

EI = ones(size(elements,1)) * E * I;

figure;
for ti = 1:length(T)
    [vph, nominal_waveno] = fn_waveguide_in_tension_dispersion(freq, E * I, T(ti), m);
    emergent_waveno = zeros(size(nominal_waveno));
    for fi = 1:length(freq)
        k = ones(size(elements,1)) * nominal_waveno(fi);
        for bi = 1:length(sleeper_nodes)
%             BC(bi) = fn_BC_values_for_lateral_mass_loading(lateral_mass, freq(fi), sleeper_nodes(bi));
            BC(bi) = fn_BC_values_for_sleeper(lateral_mass, lateral_stiffness, lateral_damping, rotational_mass, rotational_stiffness, rotational_damping, freq(fi), sleeper_nodes(bi));
        end
        [K, S] = fn_build_flex_global_matrices(nodes, elements, EI, k);
        [emergent_waveno(fi), x, mode_shape] = fn_calculate_emergent_waveno_for_model(K, S, BC, nodes, elements, k);
    end
    % figure;
    loglog(freq, 2 * pi * freq ./ emergent_waveno, [col(ti), '.-']);
    hold on;
    loglog(freq, 2 * pi * freq ./ nominal_waveno, [col(ti), '.:']);
end
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)');

figure;plot(x, [real(mode_shape), imag(mode_shape)]);