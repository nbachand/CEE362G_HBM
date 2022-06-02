function [H] = detector_H( ty, tsource, source_loc, det_loc, u, D)
%% Construct H matrix for single detector
% ty: times of measurements
% tsource: times of sources
% source_loc: source location
% det_loc: detector location
% u: wind vector,  [0,0] for example
% D: diffusion coefficient

n = length(ty);
m = length(tsource);
H = zeros(n,m); 

f = @(s_l, d_l, T) diffusion_f(s_l, d_l, T, D, u);
% f = @(s_l, d_l, T) (norm(s_l-d_l))./sqrt(4*pi*D*T.^3).*exp(-(norm(d_l-s_l-u*T).^2./(4*D*T)));

for k = 1:n
    for l = 1:m
        H(k,l) = f(source_loc, det_loc ,ty(k)-tsource(l))*1;
    end
end

