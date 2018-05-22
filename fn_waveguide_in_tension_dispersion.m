function [vph, k] = fn_waveguide_in_tension_dispersion(freq, EI, T, m)
%freq  = frequency in Hz
%EI = young's modulus * 2nd moment of area
%T = tension in Newtons
%m = mass per metre in kg


%this just implements eq. 2 in F. Chen, P.D. Wilcox / Ultrasonics 47 (2007) 111–122
w = 2*pi*freq;
vph = w .* sqrt(2*EI ./ (sqrt(T .^ 2 + 4 * m * EI * w .^ 2) - T));
k = w ./ vph;
end