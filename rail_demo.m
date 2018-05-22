clear;
close all;

%rail properties
E = 210e9;
I = 6e-6;
m = 67.7;

%tension
T = 0;

%the excitation frequency to actually use
freq = 10;

%geometry
length_of_rail = 100;
position_of_shaker = 6;
length_option = 'ends at midspans';

%sleeper details
sleeper_spacing = 0.6;
lateral_mass = 50;
lateral_stiffness = 2e5;
lateral_damping = 0.000002;
rotational_mass = 0;
rotational_stiffness = 0;
rotational_damping = 0;

%reflection coefficient at ends (0 = infinite rail, 1 = pinned support)
RC = 0;

%---------------------------------------------------------------------------
%get phase velocity and wavenumber in free rail
[true_vel_at_freq, waveno] = fn_waveguide_in_tension_dispersion(freq, E*I, T, m);
wavelength = 2 * pi / waveno;

%work out node positions
min_nodes_per_wavelength = 1; %NB you don't need to change this as it doesn't alter accuracy - only reason why elements can't be too big is because of numerical instability due to exp(kx) terms
% max_node_spacing = wavelength / min_nodes_per_wavelength;
max_node_spacing = 0.1;
[nodes, elements, sleeper_nodes, forcing_node] = fn_create_rail_mesh(length_of_rail, sleeper_spacing, position_of_shaker, max_node_spacing, length_option);

EI = ones(size(elements, 1), 1) * 1; %actual value of bending stiffness doesn't matter unless you want actual forces and moments
k = ones(size(elements, 1), 1) * waveno;

left_node = 1;
left_element = 1;
right_node = length(nodes);
right_element = size(elements, 1);

%Add sleeper BCs
for ii = 1:length(sleeper_nodes)
    BC(ii) = fn_BC_values_for_sleeper(lateral_mass, lateral_stiffness, lateral_damping, rotational_mass, rotational_stiffness, rotational_damping, freq, sleeper_nodes(ii));
end

%Add forcing BC
ii = length(BC) + 1;
BC(ii).node = forcing_node;
BC(ii).type = 'lateral forcing';
BC(ii).value = 1;

%Add end BCs
ii = length(BC) + 1;
BC(ii) = fn_BC_values_for_reflection_coefficient(RC, EI(left_element), k(left_element), left_node);
ii = length(BC) + 1;
BC(ii) = fn_BC_values_for_reflection_coefficient(RC, EI(right_element), k(right_element), right_node);


%---------------------------------------------------------------------------
%Build global matrices
[K, S] = fn_build_flex_global_matrices(nodes, elements, EI, k);

%Solve it
[u, f] = fn_flex_solver(K, BC);

%Calculate displaced shape and plot results
% xx = linspace(min(nodes), max(nodes), 500)';
pts_per_element = 10;
[xx, uu] = fn_get_displaced_shape(nodes, elements, u, S, k, pts_per_element);

figure;
% subplot(3,1,1);
plot(xx, real(uu), 'b');
hold on;
plot(xx, imag(uu), 'r');
plot(xx, abs(uu), 'k:');
plot(nodes, real(u(fn_nd_to_global_index(1:length(nodes), 1))), 'b.');
plot(nodes, imag(u(fn_nd_to_global_index(1:length(nodes), 1))), 'r.');
legend('Real', 'Imag');
xlabel('x (m)');
ylabel('Displacement');

% subplot(3,1,2);
% plot(xx, angle(uu) * 180 / pi, 'b');
% xlabel('x (m)');
% ylabel('Phase (deg)');
% 
% subplot(3,1,3);
% plot(xx, abs(2 * pi * freq ./ gradient(unwrap(angle(uu)), xx(2) - xx(1))), 'b');
% hold on;
% plot(xx, ones(size(xx)) * true_vel_at_freq, 'r:');
% ylim([0,1] * true_vel_at_freq * 1.5);
% legend('Measured', 'True');
% xlabel('x (m)');
% ylabel({'Estimated', 'velocity (m/s)'});
