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

%reflection coefficient at ends (0 = infinite rail, 1 = pinned support)
RC = 0.1;

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

left_node = 1;
forcing_node = n1;
right_node = length(nodes);

elements = [1: length(nodes) - 1; 2: length(nodes)]';
EI = ones(size(elements, 1), 1) * 1; %actual value of bending stiffness doesn't matter unless you want actual forces and moments
k = ones(size(elements, 1), 1) * waveno;

%BCs and forcing
BC(1).node = left_node;
BC(1).type = 'general damper';
BC(1).value = fn_damper_values_for_reflection_coefficient(RC, EI(1), k(1));

BC(2).node = right_node;
BC(2).type = 'general damper';
BC(2).value = fn_damper_values_for_reflection_coefficient(RC, EI(2), k(2));

BC(3).node = forcing_node;
BC(3).type = 'force';
BC(3).value = 1;

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
