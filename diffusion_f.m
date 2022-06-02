function f = diffusion_f(s_l, d_l, T, D, u)
    maskT = (T>1e-6).*T + (T<=1e-6).*1e-6;
    f = (norm(s_l-d_l))./sqrt(4*pi*D*maskT.^3).*exp(-(norm(d_l-s_l).^2./(4*D*maskT)));
    f = (T>1e-6).*f;
end

