function [Ke, Se] = fn_flex_element_matrices(x, EI, k)
%Create element stiffness matrix, Ke, for flexural element relating force
%vector, f = [Q1 Q2 M1 M2]^T, to displacement vector, 
%u = [y1 y2 theta1 theta2]^T, by f = Ke * u. 
%Matrix Se calculates shape function vector, w = [A B C D]^T, from
%w = Se * u

Te = zeros(4);
x0 = mean(x); %centre of element
% x = x - x0;
for ii = 1:2
    Te(2 * (ii - 1) + 1, :) =     fn_flex_shape_function(k, x(ii) - x0);
    Te(2 * (ii - 1) + 2, :) = k * fn_flex_shape_function(k, x(ii) - x0) .* [1i, -1i, 1, -1];
end
Se = inv(Te);
Te = zeros(4);
s = [1, -1];
for ii = 1:2
    Te(2 * (ii - 1) + 1, :) =  EI * k ^ 3 * fn_flex_shape_function(k, x(ii) - x0) .* [-1i, 1i, 1, -1] * s(ii);
    Te(2 * (ii - 1) + 2, :) = -EI * k ^ 2 * fn_flex_shape_function(k, x(ii) - x0) .* [ -1, -1, 1,  1] * s(ii);
end
Ke = Te * Se;
end