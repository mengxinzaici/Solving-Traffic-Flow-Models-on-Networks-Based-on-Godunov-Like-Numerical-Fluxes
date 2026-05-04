% 定义道路结构体
function roads = create_road(L, N, init_rho)
    roads.L = L;
    roads.N = N;
    roads.dx = L / N;
    roads.x = (0:N-1)*roads.dx + roads.dx/2;
    roads.rho = init_rho;          % 密度向量，长度 N
    roads.flux_left = 0;           % 左边界通量
    roads.flux_right = 0;          % 右边界通量
    roads.bc_type_left = 0;        % 0:外部, 1:路口
    roads.bc_type_right = 0;
end