function [veh_final, max_mass_err, max_E2, max_E3] = RKSVmain(flux_type)
    if nargin < 1
        flux_type = 'inside'; 
    end

    vmax = 1.0; rhomax = 1.0;

    % 空间网格设置
    L = 1.0; N = 150; dx = L/N;
    n_half = floor(N/2); n_rest = N - n_half;

    % k=1 谱体积法参数 (每个网格 2 个子控制体积)
    s_gl = [-1, 0, 1]; 
    h_ratios = diff(s_gl);
    
    % 1. 计算 Legendre 到 CV 均值的转换矩阵
    % P0 = 1, P1 = x (积分分别为 x 和 0.5x^2)
    M_L2C = [1, -0.5; 1, 0.5]; 
    M_C2L = inv(M_L2C)'; % 转置用于右乘
    
    % 2. 计算 CV 均值到端点的重构矩阵 R_pts
    V_pts = [1, -1; 1, 0; 1, 1]; % 勒让德基在 -1, 0, 1 处的值
    R_pts = V_pts * M_C2L'; % 大小为 3x2

    % 初始条件设置 (N x 2 矩阵)
    init_1 = [0.5 * ones(n_half, 1); 0.5 * ones(n_rest, 1)];   
    init_2 = [0.75 * ones(n_half, 1); 0.0 * ones(n_rest, 1)];   
    init_3 = [0.25 * ones(n_half, 1); 0.0 * ones(n_rest, 1)]; 
    
    rho1 = [init_1, init_1];
    rho2 = [init_2, init_2];
    rho3 = [init_3, init_3];

    alpha = [0.75; 0.25]; 

    % 时间推进设置 (k=1 匹配 SSP-RK2)
    dt = 0.5 * (dx / vmax) / 2; 
    T_end = 4.0; Nt = ceil(T_end / dt);

    plot_times = [0, 0.25, 0.5, 1.25, 2.5, 4.0];
    snapshot_idx = 1;
    rho_snapshots = cell(length(plot_times), 3);

    mass_vec = h_ratios * (dx / 2); 
    M0 = sum(rho1 * mass_vec') + sum(rho2 * mass_vec') + sum(rho3 * mass_vec');
    max_mass_err = 0; max_E2 = 0; max_E3 = 0;

    % SSP-RK2 主循环 (带严格边界限制)
    for t_step = 1:Nt
        current_time = (t_step - 1) * dt;
        
        if snapshot_idx <= length(plot_times) && current_time >= plot_times(snapshot_idx) - 1e-5
            rho_snapshots{snapshot_idx, 1} = rho1;
            rho_snapshots{snapshot_idx, 2} = rho2;
            rho_snapshots{snapshot_idx, 3} = rho3;
            snapshot_idx = snapshot_idx + 1;
        end

        % --- RK2 第一步 ---
        [L1_1, L2_1, L3_1, H_in, H_out] = get_RHS_k1(rho1, rho2, rho3, flux_type, alpha, vmax, rhomax, dx, R_pts, h_ratios);
        max_E2 = max(max_E2, abs(H_out(1) - alpha(1) * H_in(1)));
        max_E3 = max(max_E3, abs(H_out(2) - alpha(2) * H_in(1)));
        
        rho1_s1 = apply_minmod_limiter_k1(rho1 + dt * L1_1, M_C2L, M_L2C, rhomax);
        rho2_s1 = apply_minmod_limiter_k1(rho2 + dt * L2_1, M_C2L, M_L2C, rhomax);
        rho3_s1 = apply_minmod_limiter_k1(rho3 + dt * L3_1, M_C2L, M_L2C, rhomax);

        % --- RK2 第二步 ---
        [L1_2, L2_2, L3_2, H_in, H_out] = get_RHS_k1(rho1_s1, rho2_s1, rho3_s1, flux_type, alpha, vmax, rhomax, dx, R_pts, h_ratios);
        max_E2 = max(max_E2, abs(H_out(1) - alpha(1) * H_in(1)));
        max_E3 = max(max_E3, abs(H_out(2) - alpha(2) * H_in(1)));
        
        rho1 = apply_minmod_limiter_k1(0.5 * rho1 + 0.5 * (rho1_s1 + dt * L1_2), M_C2L, M_L2C, rhomax);
        rho2 = apply_minmod_limiter_k1(0.5 * rho2 + 0.5 * (rho2_s1 + dt * L2_2), M_C2L, M_L2C, rhomax);
        rho3 = apply_minmod_limiter_k1(0.5 * rho3 + 0.5 * (rho3_s1 + dt * L3_2), M_C2L, M_L2C, rhomax);
    end

    if snapshot_idx <= length(plot_times)
        rho_snapshots{snapshot_idx, 1} = rho1;
        rho_snapshots{snapshot_idx, 2} = rho2;
        rho_snapshots{snapshot_idx, 3} = rho3;
    end

    current_mass = sum(rho1 * mass_vec') + sum(rho2 * mass_vec') + sum(rho3 * mass_vec');
    max_mass_err = max(max_mass_err, abs(current_mass - M0));
    veh_final = [sum(rho1 * mass_vec'), sum(rho2 * mass_vec'), sum(rho3 * mass_vec')];
% 
%     =========================================================
%  
%     =========================================================
    X_base = zeros(4*N, 1);
    for i = 1:N
        xc = (i - 0.5) * dx;
        x_edges = xc + s_gl * (dx / 2);
        idx = (i-1)*4;
        X_base(idx+1:idx+2) = [x_edges(1); x_edges(2)]; 
        X_base(idx+3:idx+4) = [x_edges(2); x_edges(3)]; 
    end
    
    figure('Color','w','Position',[100, 100, 700, 900]);
    for k_plot = 1:length(plot_times)
        subplot(6, 1, k_plot); hold on; box on;
        tmp1 = rho_snapshots{k_plot, 1}; tmp2 = rho_snapshots{k_plot, 2}; tmp3 = rho_snapshots{k_plot, 3};
        Y1 = zeros(4*N, 1); Y2 = zeros(4*N, 1); Y3 = zeros(4*N, 1);
        flat1 = reshape(tmp1', [], 1); flat2 = reshape(tmp2', [], 1); flat3 = reshape(tmp3', [], 1);
        Y1(1:2:end) = flat1; Y1(2:2:end) = flat1;
        Y2(1:2:end) = flat2; Y2(2:2:end) = flat2;
        Y3(1:2:end) = flat3; Y3(2:2:end) = flat3;
        plot(X_base, Y1, 'r-', 'LineWidth', 1.5);
        plot(X_base + L, Y2, 'g-', 'LineWidth', 1.5);
        plot(X_base + L, Y3, 'b-', 'LineWidth', 1.5);
        xlim([0, 2]); ylim([0, 1.1]); ylabel('\rho'); title(['t = ', num2str(plot_times(k_plot))]);
        if k_plot == 1, sgtitle(['RKSV (k=1) \alpha-', flux_type, ' 通量']); legend('入路 1', '出路 2', '出路 3'); end
    end
end

%% 内部函数区：RHS 计算 (k=1)
function [L1, L2, L3, H_in, H_out] = get_RHS_k1(rho1, rho2, rho3, flux_type, alpha, vmax, rhomax, dx, R_pts, h_ratios)
    v1 = rho1 * R_pts'; v2 = rho2 * R_pts'; v3 = rho3 * R_pts'; 
    
    rho1_R_safe = max(0, min(rhomax, v1(end, 3))); 
    rho2_L_safe = max(0, min(rhomax, v2(1, 1)));   
    rho3_L_safe = max(0, min(rhomax, v3(1, 1)));   

    switch flux_type
        case 'inside'
            [H_in, H_out] = junction_flux_inside(rho1_R_safe, [rho2_L_safe; rho3_L_safe], alpha, vmax, rhomax);
        case 'outside'
            [H_in, H_out] = junction_flux_outside(rho1_R_safe, [rho2_L_safe; rho3_L_safe], alpha, vmax, rhomax);
        case 'classical'
            [H_in, H_out] = junction_flux_classical(rho1_R_safe, [rho2_L_safe; rho3_L_safe], alpha, vmax, rhomax);
    end

    L1 = compute_L_sv1(rho1, 0, H_in(1), dx, R_pts, h_ratios, vmax, rhomax);
    L2 = compute_L_sv1(rho2, H_out(1), 0, dx, R_pts, h_ratios, vmax, rhomax);
    L3 = compute_L_sv1(rho3, H_out(2), 0, dx, R_pts, h_ratios, vmax, rhomax);
end

function L_val = compute_L_sv1(rho_cv, F_left_node, F_right_node, dx, R_pts, h_ratios, vmax, rhomax)
    N = size(rho_cv, 1);
    v_pts = rho_cv * R_pts'; 
    v_pts_safe = max(0, min(rhomax, v_pts)); 

    F = zeros(N, 3);
    % 内部 CV 交界面
    F(:, 2) = greenshields_flux(v_pts_safe(:, 2), vmax, rhomax);

    % 网格外部交界面
    F(1, 1) = F_left_node;
    for i = 1:N-1
        F(i+1, 1) = godunov_flux(v_pts_safe(i, 3), v_pts_safe(i+1, 1), vmax, rhomax);
        F(i, 3) = F(i+1, 1); 
    end
    F(N, 3) = F_right_node;

    L_val = zeros(N, 2);
    for j = 1:2
        L_val(:, j) = -(2/dx) * (F(:, j+1) - F(:, j)) / h_ratios(j);
    end
end

%% 斜率限制器 (严格 TVD + 保界)
function rho_lim = apply_minmod_limiter_k1(rho, M_C2L, M_L2C, rhomax)
    U_coeffs = rho * M_C2L; 
    N = size(rho, 1);
    
    % 1. 内部网格限制
    for i = 2 : N-1
        Uc = U_coeffs(i, :); Ul = U_coeffs(i-1, :); Ur = U_coeffs(i+1, :);
        U_coeffs(i, 2) = minmod3(Uc(2), Ur(1) - Uc(1), Uc(1) - Ul(1));
    end
    
    % 2. 左边界网格严格保界 (防止漏出物理边界产生分配误差)
    slope_1 = minmod2(U_coeffs(1, 2), U_coeffs(2, 1) - U_coeffs(1, 1));
    max_allow_1 = min(rhomax - U_coeffs(1, 1), U_coeffs(1, 1));
    U_coeffs(1, 2) = sign(slope_1) * min(abs(slope_1), max_allow_1);
    
    % 3. 右边界网格严格保界
    slope_N = minmod2(U_coeffs(N, 2), U_coeffs(N, 1) - U_coeffs(N-1, 1));
    max_allow_N = min(rhomax - U_coeffs(N, 1), U_coeffs(N, 1));
    U_coeffs(N, 2) = sign(slope_N) * min(abs(slope_N), max_allow_N);
    
    rho_lim = U_coeffs * M_L2C'; 
end

function m = minmod3(a, b, c)
    if sign(a) == sign(b) && sign(b) == sign(c)
        m = sign(a) * min([abs(a), abs(b), abs(c)]);
    else
        m = 0;
    end
end

function m = minmod2(a, b)
    if sign(a) == sign(b)
        m = sign(a) * min(abs(a), abs(b));
    else
        m = 0;
    end
end