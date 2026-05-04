function [H_in, H_out] = junction_flux_classical(rho_in_end, rho_out_start, alpha, vmax, rhomax)
    n = length(rho_in_end);
    m = length(rho_out_start);
    H_out = zeros(m,1);
    H_in = zeros(n,1);
    
    rho_star = rhomax / 2;
    f_star = vmax * rho_star * (1 - rho_star/rhomax);
    
    % 计算每一条入边的实际总流出量 H_in
    for i = 1:n
        % 1. 计算入边 i 的需求 fin
        if rho_in_end(i) <= rho_star
            fin = vmax * rho_in_end(i) * (1 - rho_in_end(i)/rhomax);
        else
            fin = f_star;
        end
        
        min_flux = fin;
        
        % 2. 遍历所有出边，寻找严格遵守比例下的极小值限制
        for j = 1:m
            if rho_out_start(j) <= rho_star
                fout = f_star;
            else
                fout = vmax * rho_out_start(j) * (1 - rho_out_start(j)/rhomax);
            end
            
            % 不断更新极小值，寻找最大可能的总流量 (加 eps 防止除以 0)
            min_flux = min(min_flux, fout / (alpha(j,i) + eps));
        end
        
        H_in(i) = min_flux;
    end
    
    % 计算每一条出边的实际总流入量 H_out
    for j = 1:m
        sum_flux = 0;
        for i = 1:n
            sum_flux = sum_flux + alpha(j,i) * H_in(i);
        end
        H_out(j) = sum_flux;
    end
end