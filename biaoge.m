% generate_tables.m
clear; clc; close all;

% --- 1. 运行模拟 ---
% k=0 (一阶有限体积法)
[veh_fvm_c, err_M_fvm_c, E2_fvm_c, E3_fvm_c] = LWR1to2('classical');
[veh_fvm_o, err_M_fvm_o, E2_fvm_o, E3_fvm_o] = LWR1to2('outside');
[veh_fvm_i, err_M_fvm_i, E2_fvm_i, E3_fvm_i] = LWR1to2('inside');

% k=1 (二阶谱体积法)
[veh_rksv_c, err_M_rksv_c, E2_rksv_c, E3_rksv_c] = RKSVmain('classical');
[veh_rksv_o, err_M_rksv_o, E2_rksv_o, E3_rksv_o] = RKSVmain('outside');
[veh_rksv_i, err_M_rksv_i, E2_rksv_i, E3_rksv_i] = RKSVmain('inside');

% k=2 (三阶谱体积法) 
% 注意：之前编写的 SVM_k2_main 仅包含 inside 和 outside 两种容错通量
[veh_sv2_o, err_M_sv2_o, E2_sv2_o, E3_sv2_o] = SVM_k2_main('outside');
[veh_sv2_i, err_M_sv2_i, E2_sv2_i, E3_sv2_i] = SVM_k2_main('inside');

% 数据格式化函数
fmt_n = @(x) sprintf('%.6f', x); 
fmt_e = @(x) sprintf('%.2e', x);  

% --- 表 5.1 数据：三种节点通量最终车辆数 (仅针对 SV(k=0)) ---
data61 = {
    'Maximum possible', fmt_n(veh_fvm_c(1)), fmt_n(veh_fvm_c(2)), fmt_n(veh_fvm_c(3)), fmt_e(err_M_fvm_c);
    'α-Outside',   fmt_n(veh_fvm_o(1)), fmt_n(veh_fvm_o(2)), fmt_n(veh_fvm_o(3)), fmt_e(err_M_fvm_o);
    'α-Inside',    fmt_n(veh_fvm_i(1)), fmt_n(veh_fvm_i(2)), fmt_n(veh_fvm_i(3)), fmt_e(err_M_fvm_i)
};

% --- 表 5.2 数据：不同多项式阶数格式的质量守恒误差对比 ---
data62 = {
    'SV(k=0)',  '150', fmt_e(err_M_fvm_i);
    'RKSV(k=1)', '150', fmt_e(err_M_rksv_i);
    'RKSV(k=2)',  '150', fmt_e(err_M_sv2_i)
};

% --- 表 5.3 数据：不同阶数下的交通分配误差对比 ---
data63 = {
    'SV(k=0)',  'α-Outside', fmt_n(E2_fvm_o), fmt_n(E3_fvm_o);
    'SV(k=0)',  'α-Inside',  fmt_n(E2_fvm_i), fmt_n(E3_fvm_i);
    'RKSV(k=1)', 'α-Outside', fmt_n(E2_rksv_o), fmt_n(E3_rksv_o);
    'RKSV(k=1)', 'α-Inside',  fmt_n(E2_rksv_i), fmt_n(E3_rksv_i);
    'RKSV(k=2)', 'α-Outside', fmt_n(E2_sv2_o), fmt_n(E3_sv2_o);
    'RKSV(k=2)', 'α-Inside',  fmt_n(E2_sv2_i), fmt_n(E3_sv2_i)
};

% ==================== 界面与表格绘制 ====================
fig = figure('Name', '交通流模拟结果表格', 'Position', [200, 100, 750, 750], 'Color', 'w', 'MenuBar', 'none');

% 通用表格列宽设置
col_width = {140, 110, 110, 110, 140};

% --- 绘制 表 5.1 ---
uicontrol('Style', 'text', 'String', '表 5.1：三种节点通量最终车辆数 (基于 SV(k=0) 方法)', ...
    'Units', 'normalized', 'Position', [0.05, 0.90, 0.9, 0.04], ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

uitable(fig, 'Data', data61, ...
    'ColumnName', {'通量类型', 'Road 1', 'Road 2', 'Road 3', '总质量误差'}, ...
    'Units', 'normalized', 'Position', [0.05, 0.74, 0.9, 0.15], ...
    'RowName', [], 'ColumnWidth', col_width, 'FontSize', 11);

% --- 绘制 表 5.2 ---
uicontrol('Style', 'text', 'String', '表 5.2：各阶数谱体积法的质量守恒误差 (α-Inside 通量)', ...
    'Units', 'normalized', 'Position', [0.05, 0.62, 0.9, 0.04], ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

uitable(fig, 'Data', data62, ...
    'ColumnName', {'格式', '网格数', '最大质量误差'}, ...
    'Units', 'normalized', 'Position', [0.05, 0.45, 0.9, 0.16], ...
    'RowName', [], 'ColumnWidth', {140, 110, 220}, 'FontSize', 11);

% --- 绘制 表 5.3 ---
uicontrol('Style', 'text', 'String', '表 5.3：不同格式下交通分配误差演化对比', ...
    'Units', 'normalized', 'Position', [0.05, 0.35, 0.9, 0.04], ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

uitable(fig, 'Data', data63, ...
    'ColumnName', {'格式', '通量', 'max |E2(t)|', 'max |E3(t)|'}, ...
    'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.28], ...
    'RowName', [], 'ColumnWidth', {140, 140, 140, 140}, 'FontSize', 11);