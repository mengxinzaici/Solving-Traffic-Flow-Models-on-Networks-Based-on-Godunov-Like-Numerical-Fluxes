% 带α的Godunov通量（论文公式20）
function H = godunov_flux_alpha(rhoL, rhoR, alpha, vmax, rhomax)
    rho_star = rhomax / 2;
    f_star = vmax * rho_star * (1 - rho_star/rhomax);
    
    if rhoL <= rho_star
        fin = vmax * rhoL * (1 - rhoL/rhomax);
    else
        fin = f_star;
    end
    
    if rhoR <= rho_star
        fout = f_star;
    else
        fout = vmax * rhoR * (1 - rhoR/rhomax);
    end
    
    H = min(alpha * fin, fout);
end