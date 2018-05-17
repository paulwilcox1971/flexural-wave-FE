clear;
close all;
% clc;
%test flex stiff matrix and shape function

x = [1;8]+5;
k = 2 * pi / 7 / 3;
EI = 1;

%imposed displacements
u = [0;-i;0;1];

[Ke, Se] = fn_flex_element_matrices(x, EI, k);
Ke

return
xx = linspace(x(1), x(2), 100);
yy = fn_flex_shape_function(k, xx) * (Se * u);

figure;plot(xx, [real(yy), imag(yy)]);