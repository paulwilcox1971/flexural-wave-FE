clear;
close all;
% clc;
%test flex stiff matrix and shape function
d = 5;
nodes = [0;d/4;d];
elements = [1,2;2,3];
k = ones(size(elements,1)) * 1;
EI = ones(size(elements,1)) * 1;

[K, S] = fn_build_flex_global_matrices(nodes, elements, EI, k);

%FOR REAL ONE, ANALYSIS HAS TO BE PERFORMED ON MODIFIED K AFTER OTHER BCs
%ARE APPLIED!

%get compliance matrix
C = inv(K);

%extract first and last pairs of rows/cols
C = C([1, 2, end - 1, end], [1, 2, end - 1,end]);

C11 = C(1:2, 1:2);
C12 = C(1:2, 3:4);
C21 = C(3:4, 1:2);
C22 = C(3:4, 3:4);

%Convert to transfer matrix
invC12 = inv(C(1:2, 3:4));
T = [C22 * invC12, C21 - C22 * invC12 * C11; invC12, -invC12 * C11];T(:, 3:4) = -T(:, 3:4); %because sign convention for forces is opposite for a transfer matrix to a stiffness matrix
[V,D] = eig(T);
theta = angle(diag(D));
[dummy, n] = max(theta);
% k_extracted = theta(n) /(nodes(2) - nodes(1));
%now need to work out how many 2pi phase shifts occured over model length
left_node = 1;
right_node = length(nodes);
%BCs and forcing
BC(1).node = left_node;
BC(1).type = 'general displacement';
BC(1).value = V(1:2,n);

BC(2).node = right_node;
BC(2).type = 'general displacement';
BC(2).value = V(1:2,n) * exp(1i * theta(n));


%Solve it
[u, f] = fn_flex_solver(K, BC);

%Calculate displaced shape and plot results
xx = linspace(min(nodes), max(nodes), 500)';
uu = fn_get_displaced_shape(xx, nodes, u, S, elements, k);
theta_true = unwrap(angle(uu));
k_extracted = (theta_true(end) - theta_true(1)) / (nodes(end) - nodes(1))

figure;plot(xx, [real(uu), imag(uu)]);hold on; plot(nodes, real(u(1:2:end)), 'bo');