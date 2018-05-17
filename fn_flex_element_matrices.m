function [Ke, Se] = fn_flex_element_matrices(x, EI, k)
%Create element stiffness matrix, Ke, for flexural element relating force
%vector, f = [Q1 Q2 M1 M2]^T, to displacement vector, 
%u = [y1 y2 theta1 theta2]^T, by f = Ke * u. 
%Matrix Se calculates shape function vector, w = [A B C D]^T, from
%w = Se * u
x1 = x(1);
x2 = x(2);
Te = zeros(4);
for ii = 1:2
    Te(ii,:)   = [         exp(1i * k * x(ii)),           exp(-1i * k * x(ii)),     exp(k * x(ii)),      exp(-k * x(ii))];
    Te(ii+2,:) = [1i * k * exp(1i * k * x(ii)), -1i * k * exp(-1i * k * x(ii)), k * exp(k * x(ii)), -k * exp(-k * x(ii))];
end
Se = inv(Te);
Te = zeros(4);
for ii = 1:2
    Te(ii,:)   = [-1i * k * exp(1i * k * x(ii)), 1i * k * exp(-1i * k * x(ii)), k * exp(k * x(ii)), -k * exp(-k * x(ii))];
    Te(ii+2,:) = [         -exp(1i * k * x(ii)),         -exp(-1i * k * x(ii)),     exp(k * x(ii)),      exp(-k * x(ii))];
end
Te = -EI * k^2 * Te;
Ke = Te * Se;
end