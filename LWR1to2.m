function [veh_final, max_mass_err, max_E2, max_E3] = LWR1to2(flux_type)
% flux_type 可选参数: 'classical', 'outside', 'inside'

% clear; clc; close all;

vmax = 1.0;
rhomax = 1.0;

% 创建三条道路（入边1，出边2，出边3）
L = 1.0; N = 150; dx = L/N;
init_rho1 = [0.5*ones(N/2,1); 0.5*ones(N/2,1)];     % 入边初始密度 0.2
init_rho2 = [0.75*ones(N/2,1); 0.0*ones(N/2,1)];      % 出边2初始密度 0.2
init_rho3 = [0.25*ones(N/2,1); 0.0*ones(N/2,1)];      % 出边3初始密度 0.2

roads(1) = create_road(L, N, init_rho1);
roads(2) = create_road(L, N, init_rho2);
roads(3) = create_road(L, N, init_rho3);

% 设置边界类型
roads(1).bc_type_left = 0;   % 外部边界（入口）
roads(1).bc_type_right = 1;  % 连接到路口
roads(2).bc_type_left = 1;   % 连接到路口
roads(2).bc_type_right = 0;  % 外部边界（出口）
roads(3).bc_type_left = 1;
roads(3).bc_type_right = 0;

% 外部边界条件
rho_in = 0.0;

% 路口信息
junction.in_roads = [1];           % 入边索引
junction.out_roads = [2,3];        % 出边索引
junction.alpha = [0.75; 0.25];     % 分布系数（从入边1到出边2和3）

% 时间设置
CFL = 0.5;
dx_min = min([roads(1).dx, roads(2).dx, roads(3).dx]);   % 全局最小空间步长
dt = CFL * dx_min / vmax;
Tf = 4.0;
Nt = ceil(Tf / dt);
dt = Tf / Nt;   % 调整时间步长使总步数为整数

total_veh1 = zeros(Nt+1,1);  % Road1
total_veh2 = zeros(Nt+1,1);  % Road2
total_veh3 = zeros(Nt+1,1);  % Road3

total_veh1(1) = sum(roads(1).rho)*roads(1).dx;
total_veh2(1) = sum(roads(2).rho)*roads(2).dx;
total_veh3(1) = sum(roads(3).rho)*roads(3).dx;





% 存储总车辆数
total_veh = zeros(Nt+1,1);
total_veh(1) = sum(roads(1).rho)*roads(1).dx + sum(roads(2).rho)*roads(2).dx + sum(roads(3).rho)*roads(3).dx;



output_times = [0, 0.25, 0.5, 1.25, 2.5, 4.0];
rho_snapshots = cell(length(output_times), 3);
rho_snapshots(1,:) = {roads(1).rho, roads(2).rho, roads(3).rho};
next_out = 2;

% 记录初始全网总质量 M0
M0 = sum(roads(1).rho)*roads(1).dx + sum(roads(2).rho)*roads(2).dx + sum(roads(3).rho)*roads(3).dx;

% 初始化最大误差变量
max_mass_err = 0;
max_E2 = 0;
max_E3 = 0;





% 时间循环
for n = 1:Nt
    %收集路口边界密度
    rho_in_end = zeros(length(junction.in_roads),1);
    for k = 1:length(junction.in_roads)
        r = junction.in_roads(k);
        rho_in_end(k) = roads(r).rho(end);
    end
    rho_out_start = zeros(length(junction.out_roads),1);
    for k = 1:length(junction.out_roads)
        r = junction.out_roads(k);
        rho_out_start(k) = roads(r).rho(1);
    end
% 根据传入的通量类型选择调用的函数
    switch flux_type
        case 'inside'
            [H_in, H_out] = junction_flux_inside(rho_in_end, rho_out_start, junction.alpha, vmax, rhomax);
        case 'outside'
            [H_in, H_out] = junction_flux_outside(rho_in_end, rho_out_start, junction.alpha, vmax, rhomax);
        case 'classical'
            [H_in, H_out] = junction_flux_classical(rho_in_end, rho_out_start, junction.alpha, vmax, rhomax);
    end
    
    % 计算这一步的交通分配误差 Ej = Hj_actual - alpha_j * H_in
    E2 = H_out(1) - junction.alpha(1) * H_in(1);
    E3 = H_out(2) - junction.alpha(2) * H_in(1);
    max_E2 = max(max_E2, abs(E2));
    max_E3 = max(max_E3, abs(E3));
    % 将路口通量赋给道路的边界通量
    for i = 1:length(junction.in_roads)
        r = junction.in_roads(i);
        roads(r).flux_right = H_in(i);
    end
    for j = 1:length(junction.out_roads)
        r = junction.out_roads(j);
        roads(r).flux_left = H_out(j);
    end
    
    roads(1).flux_left = 0; 
    roads(2).flux_right = 0;
    roads(3).flux_right = 0;
    %处理外部边界通量
%     % 入边左边界
%     roads(1).flux_left = godunov_flux(rho_in, roads(1).rho(1), vmax, rhomax);
%     % 出边右边界
%     roads(2).flux_right = godunov_flux(roads(2).rho(end), roads(2).rho(end), vmax, rhomax);
%     roads(3).flux_right = godunov_flux(roads(3).rho(end), roads(3).rho(end), vmax, rhomax);
    
    %更新所有道路
    roads(1).rho = update_road(roads(1), dt, vmax, rhomax);
    roads(2).rho = update_road(roads(2), dt, vmax, rhomax);
    roads(3).rho = update_road(roads(3), dt, vmax, rhomax);
    
    total_veh1(n+1) = sum(roads(1).rho)*roads(1).dx;
    total_veh2(n+1) = sum(roads(2).rho)*roads(2).dx;
    total_veh3(n+1) = sum(roads(3).rho)*roads(3).dx;


    % 保存总车辆数
    total_veh(n+1) = sum(roads(1).rho)*roads(1).dx + sum(roads(2).rho)*roads(2).dx + sum(roads(3).rho)*roads(3).dx;


if next_out <= length(output_times) && n*dt >= output_times(next_out) - 1e-10
    rho_snapshots(next_out,:) = {roads(1).rho, roads(2).rho, roads(3).rho};
    next_out = next_out + 1;
end



% 计算当前时刻全网质量与初始质量的绝对误差
    current_mass = sum(roads(1).rho)*roads(1).dx + sum(roads(2).rho)*roads(2).dx + sum(roads(3).rho)*roads(3).dx;
    max_mass_err = max(max_mass_err, abs(current_mass - M0));



end


% figure('Color','w','Position',[100,100,400,200]);
% hold on; grid on;
% 
% x_local = linspace(dx/2, L-dx/2, N)';  
% 
% plot(x_local, roads(1).rho, 'r-', 'LineWidth',2);
% 
% plot(x_local + 1, roads(2).rho, 'g-', 'LineWidth',2);
% 
% plot(x_local + 1, roads(3).rho, 'b-', 'LineWidth',2);
% 
% xlim([0, 2]);
% ylim([0, 1.0]);
% xlabel('x');
% ylabel('\rho');
% title('入路 0-1，出路 1-2');
% %legend('Road 1','Road 2','Road 3','Location','best');
% set(gca,'FontSize',12);
% hold off;
% 
% 
% figure('Color','w');
% for k = 1:length(output_times)
%     subplot(6,1,k);
%     plot(x_local, rho_snapshots{k,1}, 'r-'); hold on;
%     plot(x_local+1, rho_snapshots{k,2}, 'g-');
%     plot(x_local+1, rho_snapshots{k,3}, 'b-');
%     xlim([0,2]); ylim([0,1]); title(['t = ', num2str(output_times(k))]);
% %     if k==1, legend('Road1','Road2','Road3'); end
%     hold off;
% end
% sgtitle('α-inside通量');

% 
% 
% figure('Color','w','Position',[100,300,500,300]);
% hold on; grid on;
% t_vec = (0:Nt)*dt;
% plot(t_vec, total_veh1, 'd-', 'LineWidth',2.5);  % Road3
% plot(t_vec, total_veh2, 'g-', 'LineWidth',2.5);  % Road2
% plot(t_vec, total_veh3, 'b-', 'LineWidth',2.5);  % Road3
% 
% xlabel('时间 t');
% ylabel('车辆总数');
% title('Road2 与 Road3 车辆数时间序列');
% legend('Road 2','Road 3','Location','best');
% ylim([0, 1]);
% set(gca,'FontSize',12);
% hold off;

veh_final = [sum(roads(1).rho)*roads(1).dx, sum(roads(2).rho)*roads(2).dx, sum(roads(3).rho)*roads(3).dx];
end % 结束 LWR1to2 函数