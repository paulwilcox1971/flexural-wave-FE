clear;
close all;

%velocity at excitation frequency will be calculated from pair of known 
%values(freq0 and vel_at_freq0) on dispersion curve assuming vel 
%proportional to sqrt(freq)
freq0 = 100;
vel_at_freq0 = 100;

%the excitation frequency to actually use
freq = 21;

%geometry
length_of_rail = 20;
position_of_shaker = 6;
position_of_mass = 12;

%reflection coefficient at ends (0 = infinite rail, 1 = pinned support)
RC = 0;

%---------------------------------------------------------------------------
%convert above into description for FE analysis
true_vel_at_freq = sqrt(freq / freq0) * vel_at_freq0;
waveno = 2 * pi * freq / true_vel_at_freq;
wavelength = 2 * pi / waveno;

%work out node positions
min_nodes_per_wavelength = 1; %NB you don't need to change this as it doesn't alter accuracy - only reason why elements can't be too big is because of numerical instability due to exp(kx) terms
n1 = max([2, ceil(position_of_shaker / wavelength * min_nodes_per_wavelength)]);
n2 = max([2, ceil((length_of_rail - position_of_shaker) / wavelength * min_nodes_per_wavelength)]);
nodes = linspace(0, position_of_shaker, n1);
nodes = [nodes(1:end - 1), linspace(position_of_shaker, length_of_rail, n2)]'; 

elements = [1: length(nodes) - 1; 2: length(nodes)]';
EI = ones(size(elements, 1), 1) * 1; %actual value of bending stiffness doesn't matter unless you want actual forces and moments
k = ones(size(elements, 1), 1) * waveno;

left_node = 1;
left_element = 1;
forcing_node = n1;
right_node = length(nodes);
right_element = size(elements, 1);


%BCs and forcing
BC(1) = fn_BC_values_for_reflection_coefficient(RC, EI(left_element), k(left_element), left_node);

BC(2) = fn_BC_values_for_reflection_coefficient(RC, EI(right_element), k(right_element), right_node);

BC(3).node = forcing_node;
BC(3).type = 'lateral forcing';
BC(3).value = 1;

%wang a mass on somewhere to see what happens
[~,ii] = min(abs(position_of_mass - nodes));
BC(4) = fn_BC_values_for_mass_loading(0, freq, ii); 

%---------------------------------------------------------------------------
%Build global matrices
[K, S] = fn_build_flex_global_matrices(nodes, elements, EI, k);

%Solve it
[u, f] = fn_flex_solver(K, BC);

%Calculate displaced shape and plot results
xx = linspace(min(nodes), max(nodes), 500)';
uu = fn_get_displaced_shape(xx, nodes, u, S, elements, k);

figure;
subplot(3,1,1);
plot(xx, real(uu), 'b');
hold on;
plot(xx, imag(uu), 'r');
plot(nodes, real(u(fn_nd_to_global_index(1:length(nodes), 1))), 'b.');
plot(nodes, imag(u(fn_nd_to_global_index(1:length(nodes), 1))), 'r.');
legend('Real', 'Imag');
xlabel('x (m)');
ylabel('Displacement');

subplot(3,1,2);
plot(xx, angle(uu) * 180 / pi, 'b');
xlabel('x (m)');
ylabel('Phase (deg)');

subplot(3,1,3);
plot(xx, abs(2 * pi * freq ./ gradient(unwrap(angle(uu)), xx(2) - xx(1))), 'b');
hold on;
plot(xx, ones(size(xx)) * true_vel_at_freq, 'r:');
ylim([0,1] * true_vel_at_freq * 1.5);
legend('Measured', 'True');
xlabel('x (m)');
ylabel({'Estimated', 'velocity (m/s)'});
