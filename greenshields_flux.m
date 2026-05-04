function flux = greenshields_flux(rho, vmax, rhomax)
    % Greenshields 流量函数 f(rho) = vmax * rho * (1 - rho/rhomax)
    flux = vmax * rho .* (1 - rho/rhomax);
end