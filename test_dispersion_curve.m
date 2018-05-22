clear;
close all;

%see if you can recover wavenumber from model!

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

%---------------------------------------------------------------------------
%Build global matrices
[K, S] = fn_build_flex_global_matrices(nodes, elements, EI, k);

