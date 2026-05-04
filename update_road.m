function rho_new = update_road(road, dt, vmax, rhomax)
    N = road.N;
    dx = road.dx;
    rho = road.rho;
    
    % 构建通量向量 H(1:N+1)
    H = zeros(N+1,1);
    H(1) = road.flux_left;
    for i = 1:N-1
        H(i+1) = godunov_flux(rho(i), rho(i+1), vmax, rhomax);
    end
    H(N+1) = road.flux_right;
    
    % 更新密度
    rho_new = rho - (dt/dx) * (H(2:end) - H(1:end-1));
    rho_new = max(0, min(rhomax, rho_new));
end