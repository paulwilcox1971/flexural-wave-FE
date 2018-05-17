function [f_applied, u_applied, k_applied] = fn_flex_BCs(BC, sz)
%creates appropriate vectors for specified BCs and forcing
dof_per_nd = 2;
f_applied = zeros(sz, 1); 
u_applied = NaN(sz, 1); 
k_applied = NaN(sz, 1); 

for ii = 1:length(BC)
    switch BC(ii).type
        case 'pinned' %no lateral displacement
            u_applied(fn_nd_to_global_index(BC(ii).node, 1)) = 0;
        case 'clamped' %no lateral displacement or rotation
            u_applied(fn_nd_to_global_index(BC(ii).node, 1)) = 0;
            u_applied(fn_nd_to_global_index(BC(ii).node, 2)) = 0;
        
        case 'lateral displacement' %imposes lateral displacement
            u_applied(fn_nd_to_global_index(BC(ii).node, 1)) = BC(ii).value;
        case 'rotaional displacement' %imposes rotation
            u_applied(fn_nd_to_global_index(BC(ii).node, 2)) = BC(ii).value;

        case 'force' %applies lateral force
            f_applied(fn_nd_to_global_index(BC(ii).node, 1)) = BC(ii).value;
        case 'moment' %applies moment
            f_applied(fn_nd_to_global_index(BC(ii).node, 2)) = BC(ii).value;
        
        case 'lateral spring' %applies lateral spring
            k_applied(fn_nd_to_global_index(BC(ii).node, 1)) = BC(ii).value;
        case 'rotational spring' %applies rotational spring
            k_applied(fn_nd_to_global_index(BC(ii).node, 2)) = BC(ii).value;
        
        case 'lateral damper' %applies lateral damping
            k_applied(fn_nd_to_global_index(BC(ii).node, 1)) = 1i * BC(ii).value;
        case 'rotational damper' %applies rotational damping
            k_applied(fn_nd_to_global_index(BC(ii).node, 2)) = 1i * BC(ii).value;
        case 'general damper' %applies lateral and rotational damping
            for jj = 1:2
                if ~isinf(BC(ii).value(jj))
                    k_applied(fn_nd_to_global_index(BC(ii).node, jj)) = 1i * BC(ii).value(jj);
                else
                    u_applied(fn_nd_to_global_index(BC(ii).node, jj)) = 0;
                    k_applied(fn_nd_to_global_index(BC(ii).node, jj)) = NaN;
                end
            end
    end
end

end