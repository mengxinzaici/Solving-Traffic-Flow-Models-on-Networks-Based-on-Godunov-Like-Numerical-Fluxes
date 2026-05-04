function H = godunov_flux(rhoL, rhoR, vmax, rhomax)
    rho_star = rhomax / 2;
    f_star = vmax * rho_star * (1 - rho_star/rhomax);
    
    % f_in
    if rhoL <= rho_star
        fin = vmax * rhoL * (1 - rhoL/rhomax);
    else
        fin = f_star;
    end
    
    % f_out
    if rhoR <= rho_star
        fout = f_star;
    else
        fout = vmax * rhoR * (1 - rhoR/rhomax);
    end
    
    H = min(fin, fout);
end