function [u, f] = fn_flex_solver(K, BC);

%BCs and forcing
[f_applied, u_applied, k_applied] = fn_flex_BCs(BC, size(K, 1));

%deal with spring BCs by adjusting K matrix and setting relevant applied
%forces to zero and displacements to NaN - these should then be treated
%same as all other nodes with specified applied forces.
kk = find(~isnan(k_applied));
K(kk, kk) = K(kk, kk) - diag(k_applied(kk)); 
f_applied(kk) = 0;
u_applied(kk) = NaN;

ii = find(~isnan(u_applied));
jj = find( isnan(u_applied));



u = zeros(size(K,1), 1);
f = zeros(size(K,1), 1);



u(jj) = K(jj, jj) \ (f_applied(jj) - K(jj, ii) * u_applied(ii));
f(ii) = K(ii, jj) * u(jj) + K(ii, ii) * u_applied(ii);
u(ii) = u_applied(ii);
end
