function [veh_final,max_mass_err,max_E2,max_E3] = SVM_k2_main(flux_type)
    if nargin<1
        flux_type='inside'; 
    end
    vmax=1.0;rhomax=1.0;
    L=1.0;N=150;dx=L/N;
    n_half=floor(N/2);n_rest=N-n_half;
    k=2;
    s_gl=[-1,-1/sqrt(5), 1/sqrt(5), 1]; 
    h_ratios = diff(s_gl); 
    M_avg = zeros(3, 3);
    for j = 1:3
        for m = 1:3
            M_avg(j, m) = (s_gl(j+1)^m - s_gl(j)^m) / (m * h_ratios(j));
        end
    end
    V_pts = zeros(4, 3);
    for pt = 1:4
        for m = 1:3
            V_pts(pt, m) = s_gl(pt)^(m-1);
        end
    end
    R_pts = V_pts / M_avg; 

    M_L2C = zeros(3, 3);
    P_ints = @(x) [x, 0.5*x^2, 0.5*x^3 - 0.5*x]; 
    for j = 1:3
        M_L2C(j, :) = (P_ints(s_gl(j+1)) - P_ints(s_gl(j))) / h_ratios(j);
    end
    M_C2L = inv(M_L2C)'; 

    rho1 = 0.5 * ones(N, 3); 
    rho2 = [0.75 * ones(n_half, 3); 0.0 * ones(n_rest, 3)];
    rho3 = [0.25 * ones(n_half, 3); 0.0 * ones(n_rest, 3)];
    alpha = [0.75; 0.25];

    dt = 0.5 * (dx / vmax) / 3; 
    T_end = 4.0; Nt = ceil(T_end / dt);

     plot_times = [0, 0.25, 0.5, 1.25, 2.5, 4.0];
    snapshot_idx = 1;
    rho_snapshots = cell(length(plot_times), 3);

    mass_vec = h_ratios * (dx / 2); 
    M0 = sum(rho1 * mass_vec') + sum(rho2 * mass_vec') + sum(rho3 * mass_vec');
    
    max_mass_err = 0; max_E2 = 0; max_E3 = 0;

    for t_step = 1:Nt
        current_time = (t_step - 1) * dt;
        
        if snapshot_idx <= length(plot_times) && current_time >= plot_times(snapshot_idx) - 1e-5
            rho_snapshots{snapshot_idx, 1} = rho1;
            rho_snapshots{snapshot_idx, 2} = rho2;
            rho_snapshots{snapshot_idx, 3} = rho3;
            snapshot_idx = snapshot_idx + 1;
        end

        % 注意：此处新增了 rhomax 传参给限制器
        [L1_1, L2_1, L3_1, H_in, H_out] = get_RHS_k2(rho1, rho2, rho3, flux_type, alpha, vmax, rhomax, dx, R_pts, h_ratios);
        max_E2 = max(max_E2, abs(H_out(1) - alpha(1) * H_in(1)));
        max_E3 = max(max_E3, abs(H_out(2) - alpha(2) * H_in(1)));
        rho1_s1 = apply_minmod_limiter(rho1 + dt * L1_1, M_C2L, M_L2C, rhomax);
        rho2_s1 = apply_minmod_limiter(rho2 + dt * L2_1, M_C2L, M_L2C, rhomax);
        rho3_s1 = apply_minmod_limiter(rho3 + dt * L3_1, M_C2L, M_L2C, rhomax);

        [L1_2, L2_2, L3_2, H_in, H_out] = get_RHS_k2(rho1_s1, rho2_s1, rho3_s1, flux_type, alpha, vmax, rhomax, dx, R_pts, h_ratios);
        max_E2 = max(max_E2, abs(H_out(1) - alpha(1) * H_in(1)));
        max_E3 = max(max_E3, abs(H_out(2) - alpha(2) * H_in(1)));
        rho1_s2 = apply_minmod_limiter(0.75 * rho1 + 0.25 * (rho1_s1 + dt * L1_2), M_C2L, M_L2C, rhomax);
        rho2_s2 = apply_minmod_limiter(0.75 * rho2 + 0.25 * (rho2_s1 + dt * L2_2), M_C2L, M_L2C, rhomax);
        rho3_s2 = apply_minmod_limiter(0.75 * rho3 + 0.25 * (rho3_s1 + dt * L3_2), M_C2L, M_L2C, rhomax);

        [L1_3, L2_3, L3_3, H_in, H_out] = get_RHS_k2(rho1_s2, rho2_s2, rho3_s2, flux_type, alpha, vmax, rhomax, dx, R_pts, h_ratios);
        max_E2 = max(max_E2, abs(H_out(1) - alpha(1) * H_in(1)));
        max_E3 = max(max_E3, abs(H_out(2) - alpha(2) * H_in(1)));
        rho1 = apply_minmod_limiter(1/3 * rho1 + 2/3 * (rho1_s2 + dt * L1_3), M_C2L, M_L2C, rhomax);
        rho2 = apply_minmod_limiter(1/3 * rho2 + 2/3 * (rho2_s2 + dt * L2_3), M_C2L, M_L2C, rhomax);
        rho3 = apply_minmod_limiter(1/3 * rho3 + 2/3 * (rho3_s2 + dt * L3_3), M_C2L, M_L2C, rhomax);
    end

    if snapshot_idx <= length(plot_times)
        rho_snapshots{snapshot_idx, 1} = rho1;
        rho_snapshots{snapshot_idx, 2} = rho2;
        rho_snapshots{snapshot_idx, 3} = rho3;
    end

    current_mass = sum(rho1 * mass_vec') + sum(rho2 * mass_vec') + sum(rho3 * mass_vec');
    max_mass_err = max(max_mass_err, abs(current_mass - M0));
    veh_final = [sum(rho1 * mass_vec'), sum(rho2 * mass_vec'), sum(rho3 * mass_vec')];

    % 画图部分
    X_base = zeros(6*N, 1);
    for i = 1:N
        xc = (i - 0.5) * dx;
        x_edges = xc + s_gl * (dx / 2);
        idx = (i-1)*6;
        X_base(idx+1:idx+2) = [x_edges(1); x_edges(2)]; 
        X_base(idx+3:idx+4) = [x_edges(2); x_edges(3)]; 
        X_base(idx+5:idx+6) = [x_edges(3); x_edges(4)]; 
    end
    
    figure('Color','w','Position',[100, 100, 400, 200]);
    for k_plot = 1:length(plot_times)
        subplot(6, 1, k_plot);
        hold on; box on;
        
        tmp1 = rho_snapshots{k_plot, 1};
        tmp2 = rho_snapshots{k_plot, 2};
        tmp3 = rho_snapshots{k_plot, 3};

        Y1 = zeros(6*N, 1); Y2 = zeros(6*N, 1); Y3 = zeros(6*N, 1);
        
        flat1 = reshape(tmp1', [], 1);
        flat2 = reshape(tmp2', [], 1);
        flat3 = reshape(tmp3', [], 1);
        
        Y1(1:2:end) = flat1; Y1(2:2:end) = flat1;
        Y2(1:2:end) = flat2; Y2(2:2:end) = flat2;
        Y3(1:2:end) = flat3; Y3(2:2:end) = flat3;
        
        plot(X_base, Y1, 'r-', 'LineWidth', 1.5);
        plot(X_base + L, Y2, 'g-', 'LineWidth', 1.5);
        plot(X_base + L, Y3, 'b-', 'LineWidth', 1.5);
        
        xlim([0, 2]); ylim([0, 1.1]);
        ylabel('\rho');
        title(['t = ', num2str(plot_times(k_plot))]);
        if k_plot == 1
            sgtitle(['SVM (k=2) \alpha-', flux_type, ' 通量']);
%             legend('入路 1', '出路 2', '出路 3', 'Location', 'northeast'); 
        end
        set(gca, 'FontSize', 10);
    end
end

%% 内部函数区
function [L1, L2, L3, H_in, H_out] = get_RHS_k2(rho1, rho2, rho3, flux_type, alpha, vmax, rhomax, dx, R_pts, h_ratios)
    v1 = rho1 * R_pts'; 
    v2 = rho2 * R_pts'; 
    v3 = rho3 * R_pts'; 
    
    % 现在边界被严格保界了，这里的截断其实可以去掉了，但留着作为双重物理兜底也无妨
    rho1_R_safe = max(0, min(rhomax, v1(end, 4))); 
    rho2_L_safe = max(0, min(rhomax, v2(1, 1)));   
    rho3_L_safe = max(0, min(rhomax, v3(1, 1)));   

    switch flux_type
        case 'inside'
            [H_in, H_out] = junction_flux_inside(rho1_R_safe, [rho2_L_safe; rho3_L_safe], alpha, vmax, rhomax);
        case 'outside'
            [H_in, H_out] = junction_flux_outside(rho1_R_safe, [rho2_L_safe; rho3_L_safe], alpha, vmax, rhomax);
    end

    L1 = compute_L_sv2(rho1, 0, H_in(1), dx, R_pts, h_ratios, vmax, rhomax);
    L2 = compute_L_sv2(rho2, H_out(1), 0, dx, R_pts, h_ratios, vmax, rhomax);
    L3 = compute_L_sv2(rho3, H_out(2), 0, dx, R_pts, h_ratios, vmax, rhomax);
end

function L_val = compute_L_sv2(rho_cv, F_left_node, F_right_node, dx, R_pts, h_ratios, vmax, rhomax)
    N = size(rho_cv, 1);
    v_pts = rho_cv * R_pts'; 
    
    v_pts_safe = max(0, min(rhomax, v_pts)); 

    F = zeros(N, 4);
    for j = 2:3
        F(:, j) = greenshields_flux(v_pts_safe(:, j), vmax, rhomax);
    end

    F(1, 1) = F_left_node;
    for i = 1:N-1
        F(i+1, 1) = godunov_flux(v_pts_safe(i, 4), v_pts_safe(i+1, 1), vmax, rhomax);
        F(i, 4) = F(i+1, 1); 
    end
    F(N, 4) = F_right_node;

    L_val = zeros(N, 3);
    for j = 1:3
        L_val(:, j) = -(2/dx) * (F(:, j+1) - F(:, j)) / h_ratios(j);
    end
end

%% 终极形态：TVD 单侧保界斜率限制器 (Bound-Preserving Limiter)
function rho_lim = apply_minmod_limiter(rho, M_C2L, M_L2C, rhomax)
    U_coeffs = rho * M_C2L; 
    N = size(rho, 1);
    
    % 1. 内部网格：常规分层 Minmod 限制
    for i = 2 : N-1
        Uc = U_coeffs(i, :); Ul = U_coeffs(i-1, :); Ur = U_coeffs(i+1, :);
        
        a2 = Uc(3); b2 = Ur(2) - Uc(2); c2 = Uc(2) - Ul(2);
        mod_U2 = minmod3(a2, b2, c2);
        
        if mod_U2 ~= Uc(3)
            U_coeffs(i, 3) = mod_U2;
            a1 = Uc(2); b1 = Ur(1) - Uc(1); c1 = Uc(1) - Ul(1);
            U_coeffs(i, 2) = minmod3(a1, b1, c1);
        end
    end
    
    % 2. 边界网格：单侧 TVD + 严格物理保界 (Zhang-Shu 思想)
    % 不再粗暴归零，而是允许它拥有能承受的最大安全斜率，保留极端的物理梯度信息！
    
    % 左边界网格 (i=1) 
    U_coeffs(1, 3) = 0; % 仍需清除曲率防越界
    % 使用右邻居计算单侧前向斜率
    slope_1 = minmod2(U_coeffs(1, 2), U_coeffs(2, 1) - U_coeffs(1, 1));
    % 计算物理边界允许的最大斜率容量
    max_allow_1 = min(rhomax - U_coeffs(1, 1), U_coeffs(1, 1));
    % 保界截断斜率（绝不截断均值！）
    U_coeffs(1, 2) = sign(slope_1) * min(abs(slope_1), max_allow_1);
    
    % 右边界网格 (i=N) 
    U_coeffs(N, 3) = 0; 
    % 使用左邻居计算单侧后向斜率
    slope_N = minmod2(U_coeffs(N, 2), U_coeffs(N, 1) - U_coeffs(N-1, 1));
    % 计算物理边界允许的最大斜率容量
    max_allow_N = min(rhomax - U_coeffs(N, 1), U_coeffs(N, 1));
    % 保界截断斜率
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