clear;
close all;

%This scripts demonstrates use of FE code to model forced excitation of a
%rail with sleepers. Length of rail has finite length, but if end reflection
%coefficient RC = 0 it will be part of an infinite rail. 

%The code also calls a function to work out the effective wave dispersion
%(i.e. complex wavenumber) for waves in the rail with sleepers, as this is
%required to obtain the non-reflecting end conditions.

%rail properties
E = 210e9;
I = 6e-6;
m = 67.7;

%tension
T = 0;

%the excitation frequency to actually use
freq = 20;

%geometry
length_of_rail = 20;
position_of_shaker = 6;
length_option = 'ends at midspans';

%sleeper details
sleeper_spacing = 0.6;
lateral_mass = 40;
lateral_stiffness = 0;
lateral_damping = 0;
rotational_mass = 0;
rotational_stiffness = 0;
rotational_damping = 0;

%reflection coefficient at ends (0 = infinite rail, 1 = pinned support)
RC = 0;

%---------------------------------------------------------------------------
%Get phase velocity and wavenumber in free rail under tension - this
%describes the behaviour of waves in the rail between sleepers
[free_rail_vph, free_rail_waveno] = fn_waveguide_in_tension_dispersion(freq, E * I, T, m); %this is just the analytic function for dispersion curves based on the Euler Bernouilli beam equation in Chen and Wilcox paper
free_rail_wavelength = 2 * pi / free_rail_waveno;
min_nodes_per_wavelength = 1; %NB you don't need to change this as it doesn't alter accuracy - only reason why elements can't be too big is because of numerical instability due to exp(kx) terms
max_node_spacing = free_rail_wavelength / min_nodes_per_wavelength;

%---------------------------------------------------------------------------
%First do a dummy mesh with one sleeper for modal analysis to get effective 
%dispersion relation. This is needed to get non-reflecting BCs in model of 
%finite length of rail
no_sleepers = 1;
[nodes, elements, sleeper_nodes, forcing_node] = fn_create_rail_mesh(no_sleepers * sleeper_spacing, sleeper_spacing, [], max_node_spacing, 'ends at midspans');
clear('BC');
for bi = 1:length(sleeper_nodes)
    BC(bi) = fn_BC_values_for_sleeper(lateral_mass, lateral_stiffness, lateral_damping, rotational_mass, rotational_stiffness, rotational_damping, freq, sleeper_nodes(bi));
end
[K, S] = fn_build_flex_global_matrices(nodes, elements, E * I * ones(size(elements, 1), 1), free_rail_waveno * ones(size(elements, 1), 1));
[emergent_waveno, emergent_mode_shape] = fn_calculate_emergent_waveno_for_model(K, S, BC, nodes, elements, free_rail_waveno * ones(size(elements, 1)));

%---------------------------------------------------------------------------
%Now the actual model of a finite length of rail
[nodes, elements, sleeper_nodes, forcing_node] = fn_create_rail_mesh(length_of_rail, sleeper_spacing, position_of_shaker, max_node_spacing, length_option);

EI = ones(size(elements, 1), 1) * E * I; 
k = ones(size(elements, 1), 1) * free_rail_waveno;

left_node = 1;
left_element = 1;
right_node = length(nodes);
right_element = size(elements, 1);

%Add sleeper BCs
clear('BC');
for bi = 1:length(sleeper_nodes)
    BC(bi) = fn_BC_values_for_sleeper(lateral_mass, lateral_stiffness, lateral_damping, rotational_mass, rotational_stiffness, rotational_damping, freq, sleeper_nodes(bi));
end

%Add forcing BC
bi = length(BC) + 1;
BC(bi).node = forcing_node;
BC(bi).type = 'lateral forcing';
BC(bi).value = 1;

%Add end BCs
bi = length(BC) + 1;
BC(bi) = fn_BC_values_for_reflection_coefficient(RC, E * I, emergent_mode_shape, left_node);
bi = length(BC) + 1;
BC(bi) = fn_BC_values_for_reflection_coefficient(RC, E * I, emergent_mode_shape, right_node);

%Solve the forced problem
[K, S] = fn_build_flex_global_matrices(nodes, elements,  E * I * ones(size(elements, 1), 1), free_rail_waveno * ones(size(elements, 1)));
[u, f] = fn_flex_solver(K, BC);

%--------------------------------------------------------------------------
%Calculate displaced shape and plot results
pts_per_element = 10;
[xx, uu] = fn_get_displaced_shape(nodes, elements, u, S, k, pts_per_element);

figure;
h_real = plot(xx, real(uu), 'b');
hold on;
h_imag = plot(xx, imag(uu), 'r');
h_abs = plot(xx, abs(uu), 'k:');
xlim([nodes(1), nodes(end)]);
yy = ylim;
h_force = plot([1,1] * nodes(forcing_node), yy, 'g'); %mark forcing location
plot(nodes(sleeper_nodes), real(u(fn_nd_to_global_index(sleeper_nodes, 1))), 'bs');
plot(nodes(sleeper_nodes), imag(u(fn_nd_to_global_index(sleeper_nodes, 1))), 'r.');
xlabel('x (m)');
ylabel('Displacement');
legend([h_real, h_imag, h_abs, h_force], {'Real', 'Imag', 'Abs','Force'});

