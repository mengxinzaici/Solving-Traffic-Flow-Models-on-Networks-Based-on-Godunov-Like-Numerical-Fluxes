%α-inside路口通量（含α在内部）
function [H_in, H_out] = junction_flux_inside(rho_in_end, rho_out_start, alpha, vmax, rhomax)
    n = length(rho_in_end);
    m = length(rho_out_start);
    H_out = zeros(m,1);
    H_in = zeros(n,1);
    
    for j = 1:m
        sum_flux = 0;
        for i = 1:n
            sum_flux = sum_flux + godunov_flux_alpha(rho_in_end(i), rho_out_start(j), alpha(j,i), vmax, rhomax);
        end
        H_out(j) = sum_flux;
    end
    
    for i = 1:n
        sum_flux = 0;
        for j = 1:m
            sum_flux = sum_flux + godunov_flux_alpha(rho_in_end(i), rho_out_start(j), alpha(j,i), vmax, rhomax);
        end
        H_in(i) = sum_flux;
    end
end