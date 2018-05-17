clear;
close all;


nodes = [0; 4; 20];
elements = [1, 2; 2, 3];
EI = ones(size(elements, 1), 1) * 432692;
k = ones(size(elements, 1), 1) * 2 * pi /7;

%need to do BCs properly - needed constrained, forced, spring for each DOF

%---------------------------------------------------------------------------
%useful numbers
no_els = size(elements, 1);
dof_per_nd = 2;
nds_per_el = 2;
dof_per_el = nds_per_el * dof_per_nd;
const_per_el = 4;
no_nds = size(nodes, 1);

% fn_global_nd_indices = @(nd_no, dof_per_nd) [(nd_no - 1)* dof_per_nd + 1: nd_no * dof_per_nd];
% fn_global_el_indices = @(el_no, const_per_el) [(el_no - 1)* const_per_el + 1: el_no * const_per_el];

%get element matrices
Ke = zeros(dof_per_el, dof_per_el, no_els);
Se = zeros(dof_per_el, dof_per_el, no_els);
for ei = 1:no_els
    [Ke(:, :, ei), Se(:, :, ei)] = fn_flex_element_matrices([nodes(elements(ei,:))], EI(ei), k(ei));
end

%build global matrices
K = zeros(no_nds * dof_per_nd);
S = zeros(const_per_el * no_els, no_nds * dof_per_nd);
for ei = 1:no_els
    jj = [fn_global_nd_indices(elements(ei, 1)), fn_global_nd_indices(elements(ei, 2))];
    K(jj, jj) = K(jj, jj) + Ke(:, :, ei);
    ii = fn_global_el_indices(ei);
    S(ii, jj) = Se(:, :, ei);
end

%solve - bodge from here on, but working!
u_applied = NaN(size(K,1), 1);
u_applied(1) = 0; %pinned end
u_applied(5) = 0; %pinned end
f_applied = zeros(size(K,1), 1);
f_applied(3) = 1; %forcing node 2

ii = find(~isnan(u_applied));
jj = find( isnan(u_applied));

u = zeros(size(K,1), 1);
f = zeros(size(K,1), 1);
u(jj) = K(jj, jj) \ (f_applied(jj) - K(jj, ii) * u_applied(ii));
f(ii) = K(ii, jj) * u(jj) + K(ii, ii) * u_applied(ii);
u(ii) = u_applied(ii);

ABCD = S * u;

x = linspace(min(nodes), max(nodes), 500)';
[y, y_at_nodes] = fn_get_displaced_shape(x, nodes, elements, ABCD, k);

figure; 
plot(x, real(y), 'b');
hold on;
plot(x, imag(y), 'r');
plot(nodes, real(y_at_nodes), 'b.');
plot(nodes, imag(y_at_nodes), 'r.');