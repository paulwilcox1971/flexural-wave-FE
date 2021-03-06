function [emergent_waveno, emergent_mode_shape, x, mode_shape] = fn_calculate_emergent_waveno_for_model(K, S, BC, nodes, elements, k)
left_node = 1;
right_node = length(nodes);
d = abs(nodes(right_node) - nodes(left_node));

%check no BCs applied at first or last node
if any([BC(:).node] == 1) || any([BC(:).node] == length(nodes))
    error('End nodes must not have BCs for emergent wavenumber calculation');
end

%Deal with stiffness BCs by including in stiffness matrix
[f_applied, u_applied, k_applied] = fn_flex_BCs(BC, size(K, 1));

%check no BCs involved displacement or forcing
if any(f_applied) || ~all(isnan(u_applied))
    error('Cannot have forcing or displacement input BCs for emergent wavenumber calculation');
end

kk = find(~isnan(k_applied));
K_mod = K;
K_mod(kk, kk) = K_mod(kk, kk) - diag(k_applied(kk)); 

%get compliance matrix
C = inv(K_mod);

%extract first and last pairs of rows/cols
C = C([1, 2, end - 1, end], [1, 2, end - 1,end]);

C11 = C(1:2, 1:2);
C12 = C(1:2, 3:4);
C21 = C(3:4, 1:2);
C22 = C(3:4, 3:4);

%Convert to transfer matrix
invC12 = inv(C(1:2, 3:4));
T = [C22 * invC12, C21 - C22 * invC12 * C11; invC12, -invC12 * C11];
T(:, 3:4) = -T(:, 3:4); %because sign convention for forces is opposite for a transfer matrix to a stiffness matrix

%Eigenvalue analysis to get modeshapes and (wrapped) phase between ends of
%model
[V,D] = eig(T);
D = diag(D);
% theta = angle(diag(D));

%run actual model for each possible phase to figure out true wavenumber of
%propagating wave in positive direction
pts_per_element = 4;
i1 = length(BC) + 1;
i2 = length(BC) + 2;
BC(i1).node = left_node;
BC(i1).type = 'combined displacement';
BC(i2).node = right_node;
BC(i2).type = 'combined displacement';
for n = 1:4
    %BCs based on phase shift and mode shapes for nth Eigenvalue
    BC(i1).value = V(1:2,n);
    BC(i2).value = V(1:2,n) * D(n);%exp(1i * theta(n));
    
    %Solve it
    [u, f] = fn_flex_solver(K, BC);
    
    %Calculate displaced shape and plot results
    [x, mode_shape(:,n)] = fn_get_displaced_shape(nodes, elements, u, S, k, pts_per_element);
    theta_true = unwrap(angle(mode_shape(:,n)));
    theta_true = theta_true(end) - theta_true(1);
    emergent_waveno(n) = theta_true / d + 1i * log(abs(D(n))) / d;
end

%take the one with the biggest positive real wavenumber
[dummy, n] = max(real(emergent_waveno));
emergent_waveno = emergent_waveno(n);
mode_shape = mode_shape(:, n);
emergent_mode_shape = V(:, n);

end