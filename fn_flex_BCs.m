function [f_applied, u_applied, k_applied] = fn_flex_BCs(BC, sz)
%creates appropriate vectors for specified BCs and forcing
%real or complex values specified in BCs are applied exactly as they are
%stated - e.g. there is no multiplication by freq to get viscous damping.
%Any relevant calculations should be done when specifying values to use.

f_applied = zeros(sz, 1); 
u_applied = NaN(sz, 1); 
k_applied = NaN(sz, 1); 

for ii = 1:length(BC)
    switch BC(ii).type
        case 'lateral displacement' %imposes lateral displacement
            u_applied(fn_nd_to_global_index(BC(ii).node, 1)) = BC(ii).value;
        case 'rotaional displacement' %imposes rotation
            u_applied(fn_nd_to_global_index(BC(ii).node, 2)) = BC(ii).value;
            case 'combined displacement' %applies lateral displacement and rotation
            for jj = 1:2
                u_applied(fn_nd_to_global_index(BC(ii).node, jj)) = BC(ii).value(jj);
            end

        case 'lateral forcing' %applies lateral force
            f_applied(fn_nd_to_global_index(BC(ii).node, 1)) = BC(ii).value;
        case 'rotational forcing' %applies moment
            f_applied(fn_nd_to_global_index(BC(ii).node, 2)) = BC(ii).value;
        case 'general forcing' %applies lateral force and moment
            for jj = 1:2
                f_applied(fn_nd_to_global_index(BC(ii).node, jj)) = BC(ii).value(jj);
            end
        
        case 'lateral stiffness' %applies lateral spring
            k_applied(fn_nd_to_global_index(BC(ii).node, 1)) = BC(ii).value;
        case 'rotational stiffness' %applies rotational spring
            k_applied(fn_nd_to_global_index(BC(ii).node, 2)) = BC(ii).value;
        case 'combined stiffness' %applies lateral and rotational spring
            for jj = 1:2
                k_applied(fn_nd_to_global_index(BC(ii).node, jj)) = BC(ii).value(jj);
            end        
    end
    
    %deal with infinite stiffness values by instead imposing zero
    %displacement condition instead
    for jj = 1:2
        if isinf(k_applied(fn_nd_to_global_index(BC(ii).node, jj)))
            u_applied(fn_nd_to_global_index(BC(ii).node, jj)) = 0;
            k_applied(fn_nd_to_global_index(BC(ii).node, jj)) = NaN;
        end
    end
    
end

end