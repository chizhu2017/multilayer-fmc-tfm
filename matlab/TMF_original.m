
%% ================================================================
% 多层介质超声 FMC-TFM 反演成像
% 原始算法保持版 + 低内存优化
%
% 第一层：管片 350 mm
% 第二层：注浆 200 mm
% 第三层：围岩  50 mm
%
% 保留原始方法：
% 1. 原始回波最大正峰值识别
% 2. 原始表面波截断
% 3. 原始时间零点校正
% 4. 第一层直线走时及几何加权
% 5. 第二、三层离散最短走时
% 6. 原始实数相干叠加
% 7. 原始Hilbert包络及jet显示
%
% 优化：
% 1. 逐文件读取，减少FMC内存占用
% 2. 走时查找表代替重复界面搜索
% 3. 自动读取阵元数量和采样间隔
% 4. single精度计算
%
% 注意：
% 适用于水平分层、整数mm阵元坐标
% TXT文件：1.txt、2.txt、...、N.txt
% 第一列为时间(s)，其余列为回波数据
%% ================================================================

clc;
clear;
close all;

%% ================================================================
% 1. 参数设置
%% ================================================================

file_path = ...
'data';

% 原始参数
cut_point = 20;

x_loc = 50;
gap = 50;

Nx = 600;

% 三层厚度 mm
h1 = 350;
h2 = 200;
h3 = 50;

Nz = h1+h2+h3;

% 三层声速 m/s
v1 = single(4200);
v2 = single(2800);
v3 = single(4600);

% 界面位置
z_int1 = h1;
z_int2 = h1+h2;

% 时间零点微调
t0_fine = single(0);

% 成像显示
color_ratio = 0.7;

% 是否显示界面
show_interface = false;

% 是否保存数据
save_result = true;

%% ================================================================
% 2. 自动识别文件与数据
%% ================================================================

fprintf('\n========================================\n');
fprintf('    Multi-layer FMC-TFM Reconstruction\n');
fprintf('========================================\n');

files = dir(fullfile(file_path,'*.txt'));

sensor_num = numel(files);

if sensor_num==0
    error('指定文件夹中没有找到TXT文件。');
end

% 使用原始数字文件顺序
for k = 1:sensor_num

    filename = fullfile( ...
        file_path,sprintf('%d.txt',k));

    if ~isfile(filename)

        error('缺少文件：%d.txt',k);

    end

end

% 读取第一个文件
D0 = readmatrix(fullfile(file_path,'1.txt'));

Nt = size(D0,1);

nRec = size(D0,2)-1;

if nRec~=sensor_num

    error('发射文件数量与接收通道数不一致。');

end

% 与原始代码一致
time = single(D0(:,1));

if any(diff(time)<=0)
    error('时间列必须严格递增。');
end

fs = round(mean(1./diff(time)));

dt = single(1/fs);

fprintf('\n数据参数：\n');
fprintf('阵元数量：%d\n',sensor_num);
fprintf('FMC通道数量：%d\n',sensor_num^2);
fprintf('采样点数：%d\n',Nt);
fprintf('采样间隔：%.6f us\n',double(dt)*1e6);

%% ================================================================
% 3. 探头坐标
%% ================================================================

% 完全保留原始公式
probeX = single( ...
    (x_loc-gap)+(1:sensor_num)*gap);

fprintf('\n探头坐标：\n');
fprintf('起点：%.2f mm\n',probeX(1));
fprintf('终点：%.2f mm\n',probeX(end));
fprintf('间距：%.2f mm\n',gap);

% 查找表要求整数坐标
if any(abs(double(probeX)-round(double(probeX)))>1e-6)

    error('本优化版要求阵元位于整数mm坐标。');

end

if any(probeX<1) || any(probeX>Nx)

    error(['探头超出1~%d mm范围。' ...
        '请核对x_loc、gap及阵元数量。'],Nx);

end

%% ================================================================
% 4. 第一次读取：峰值索引
%
% 原始方法：
% [~,surf_arrival_index] = max(Z_raw)
%
% 这里只存储峰值位置，不保存整个FMC矩阵
%% ================================================================

fprintf('\n第一次读取：计算原始峰值...\n');

peakIndex = zeros( ...
    sensor_num,sensor_num,'single');

for k = 1:sensor_num

    if k==1

        D = D0;

    else

        D = readmatrix(fullfile( ...
            file_path,sprintf('%d.txt',k)));

    end

    if size(D,1)~=Nt || ...
       size(D,2)~=sensor_num+1

        error('文件%d.txt尺寸不一致。',k);

    end

    if max(abs(single(D(:,1))-time))> ...
            single(0.01)*dt

        error('文件%d.txt时间轴不一致。',k);

    end

    % 原始最大正峰值
    [~,idx] = max(single(D(:,2:end)),[],1);

    peakIndex(k,:) = single(idx);

end

fprintf('峰值索引计算完成。\n');

%% ================================================================
% 5. 原始时间零点校正
%% ================================================================

fprintf('\n计算原始时间零点...\n');

t0_candidate = zeros( ...
    sensor_num,sensor_num,'single');

for k = 1:sensor_num

    for j = 1:sensor_num

        % 保留原始发射接收编号关系
        tx = j;
        rx = k;

        x_t = single( ...
            (x_loc-gap)+tx*gap);

        x_r = single( ...
            (x_loc-gap)+rx*gap);

        direct_dist = abs(x_t-x_r);

        direct_time = direct_dist*1e-3/v1;

        direct_index = direct_time/dt;

        t0_candidate(k,j) = ...
            peakIndex(k,j)-single(direct_index);

    end

end

valid_t0 = isfinite(t0_candidate) & ...
           t0_candidate>0;

if ~any(valid_t0(:))

    error('无法计算有效时间零点。');

end

t0_index = median(t0_candidate(valid_t0));

t0_index = single(t0_index)+t0_fine;

fprintf('t0_index = %.6f\n',double(t0_index));

%% ================================================================
% 6. 建立原始成像网格
%% ================================================================

fprintf('\n建立成像网格...\n');

x_axis = single(1:Nx);

z1_range = single(1:z_int1-1);

z2_range = single(z_int1:z_int2-1);

z3_range = single(z_int2:Nz);

% 第一层网格
[X1,Z1] = meshgrid(x_axis,z1_range);

N1 = numel(z1_range);
N2 = numel(z2_range);
N3 = numel(z3_range);

fprintf('管片层：%d × %d\n',N1,Nx);
fprintf('注浆层：%d × %d\n',N2,Nx);
fprintf('围岩层：%d × %d\n',N3,Nx);

%% ================================================================
% 7. 第二、三层原始离散走时查找表
%
% 原始算法：
% 在界面x=1:600 mm搜索最短传播时间
%
% 优化原理：
% 对水平分层结构，走时只取决于
% 发射点与成像点的水平距离和深度。
%
% 因此先对水平整数偏移d建立走时表，
% 然后通过查表恢复每个阵元的传播时间。
%
% 不采用Snell连续近似，
% 保留原始1 mm离散搜索精度。
%% ================================================================

fprintf('\n预计算原始离散走时查找表...\n');

tic;

% 最大水平偏移
max_dx = max(abs([ ...
    double(probeX(:))-1; ...
    double(probeX(:))-Nx]));

max_dx = ceil(max_dx);

% 第二层深度差
dz2 = single(z2_range-z_int1);
dz2 = dz2(:);

% 第三层深度差
dz3 = single(z3_range-z_int2);
dz3 = dz3(:);

% 查找表
LUT2 = zeros(N2,max_dx+1,'single');

LUT3 = zeros(N3,max_dx+1,'single');

% 第一界面到第二界面的最短走时
T12 = zeros(1,max_dx+1,'single');

%% ---------- 7.1 第一至第二界面的走时 ----------

for d = 0:max_dx

    u = single(0:d);

    rem_dx = single(d)-u;

    % 第一层传播
    D1 = sqrt( ...
        u.^2+single(h1)^2);

    % 第二层传播
    D2 = sqrt( ...
        rem_dx.^2+single(h2)^2);

    % 原始离散最短路径
    T12(d+1) = min(D1/v1+D2/v2);

end

%% ---------- 7.2 第二层和第三层走时 ----------

for d = 0:max_dx

    % 整数mm界面搜索位置
    u = single(0:d);

    rem_dx = single(d)-u;

    % ---------- 第二层 ----------

    % 管片内传播路径
    D1 = sqrt( ...
        u.^2+single(h1)^2);

    T1 = D1/v1;

    % 注浆层内传播路径
    D2 = sqrt( ...
        rem_dx.^2+dz2.^2);

    % 原始离散最短走时
    T = T1+D2/v2;

    LUT2(:,d+1) = min(T,[],2);

    % ---------- 第三层 ----------

    % 前两层到第二界面的最短走时
    Tpre = T12(1:d+1);

    % 第二界面至地层成像点
    D3 = sqrt( ...
        rem_dx.^2+dz3.^2);

    T = Tpre+D3/v3;

    LUT3(:,d+1) = min(T,[],2);

end

clear T D1 D2 D3 T1 Tpre T12;

fprintf('离散走时查找表完成。\n');
fprintf('计算耗时：%.2f s\n',toc);

%% ================================================================
% 8. 为每个阵元建立传播距离和走时缓存
%
% 第一层存距离
% 第二、三层存原始离散走时
%
% 不建立大型二维界面搜索矩阵
%% ================================================================

fprintf('\n构建阵元传播缓存...\n');

n = sensor_num;

% 第一层直线距离
Dist1 = zeros(N1,Nx,n,'single');

% 第二层走时
Travel2 = zeros(N2,Nx,n,'single');

% 第三层走时
Travel3 = zeros(N3,Nx,n,'single');

for k = 1:n

    x_t = probeX(k);

    % ---------- 第一层 ----------

    Dist1(:,:,k) = sqrt( ...
        (X1-x_t).^2+Z1.^2);

    % ---------- 第二、三层 ----------

    dx = abs(double(x_axis)-double(x_t));

    % 整数mm偏移
    idx = round(dx)+1;

    Travel2(:,:,k) = LUT2(:,idx);

    Travel3(:,:,k) = LUT3(:,idx);

end

clear LUT2 LUT3;

fprintf('传播缓存构建完成。\n');

memory_MB = ...
    (numel(Dist1)+numel(Travel2)+numel(Travel3)) ...
    *4/1024^2;

fprintf('传播缓存内存约：%.2f MiB\n',memory_MB);

%% ================================================================
% 9. 原始TFM反演
%
% 保留：
% - 原始表面波截断
% - 原始时间零点
% - 第一层几何加权
% - 第二、三层相干叠加
%
% 不增加滤波、归一化和图像增强
%% ================================================================

fprintf('\n开始TFM反演...\n');

TFM_img1 = zeros(N1,Nx,'single');

TFM_img2 = zeros(N2,Nx,'single');

TFM_img3 = zeros(N3,Nx,'single');

tic;

for k = 1:n

    %% ---------- 读取当前TXT ----------

    if k==1

        D = D0;

    else

        D = readmatrix(fullfile( ...
            file_path,sprintf('%d.txt',k)));

    end

    if size(D,1)~=Nt || size(D,2)~=n+1

        error('文件%d.txt尺寸不一致。',k);

    end

    if max(abs(single(D(:,1))-time))> ...
            single(0.01)*dt

        error('文件%d.txt时间轴不一致。',k);

    end

    signals = single(D(:,2:end));

    clear D;

    %% ---------- 遍历当前TXT的接收通道 ----------

    for j = 1:n

        % 保留原始编号关系
        tx = j;
        rx = k;

        % 原始有符号回波
        sig = signals(:,j);

        %% ---------- 原始表面波截断 ----------

        surf_arrival_index = peakIndex(k,j);

        cutEnd = floor( ...
            double(surf_arrival_index)+cut_point);

        cutEnd = max(0,min(Nt,cutEnd));

        if cutEnd>0

            sig(1:cutEnd) = 0;

        end

        %% ====================================================
        % 第一层：管片层
        %% ====================================================

        dist_tx = Dist1(:,:,tx);

        dist_rx = Dist1(:,:,rx);

        % 原始总传播距离
        dist = dist_tx+dist_rx;

        % 原始传播时间
        t = dist*1e-3/v1;

        % 原始采样索引
        it = t/dt+t0_index;

        % 原始线性插值
        S1 = tfmInterp(sig,it);

        % 原始几何加权
        apod = ...
            (Z1./dist_tx).*(Z1./dist_rx);

        apod(~isfinite(apod)) = 0;

        % 相干叠加
        TFM_img1 = TFM_img1+S1.*apod;

        %% ====================================================
        % 第二层：注浆层
        %% ====================================================

        % 原始离散最短传播时间
        t_total = ...
            Travel2(:,:,tx)+Travel2(:,:,rx);

        % 原始时间索引
        it = t_total*1e-3/dt+t0_index;

        % 原始插值
        S2 = tfmInterp(sig,it);

        % 原始相干叠加
        TFM_img2 = TFM_img2+S2;

        %% ====================================================
        % 第三层：围岩层
        %% ====================================================

        t_total = ...
            Travel3(:,:,tx)+Travel3(:,:,rx);

        it = t_total*1e-3/dt+t0_index;

        S3 = tfmInterp(sig,it);

        TFM_img3 = TFM_img3+S3;

    end

    clear signals;

    %% ---------- 计算进度 ----------

    if k==1 || mod(k,5)==0 || k==n

        fprintf('进度：%d/%d，耗时%.2f s\n', ...
            k,n,toc);

    end

end

fprintf('\nTFM相干叠加完成。\n');

clear Dist1 Travel2 Travel3 D0;

%% ================================================================
% 10. 合并原始三层成像结果
%% ================================================================

TFM_img = zeros(Nz,Nx,'single');

TFM_img(1:z_int1-1,:) = TFM_img1;

TFM_img(z_int1:z_int2-1,:) = TFM_img2;

TFM_img(z_int2:Nz,:) = TFM_img3;

clear TFM_img1 TFM_img2 TFM_img3;

%% ================================================================
% 11. 原始Hilbert包络
%
% 与最初代码完全保持相同的处理方式
%% ================================================================

fprintf('计算Hilbert包络...\n');

TFM_env = abs(hilbert(TFM_img));

TFM_env = single(TFM_env);

max_value = max(TFM_env(:));

if max_value<=0

    error('成像结果无有效信号。');

end

%% ================================================================
% 12. 原始成像显示
%
% 保留：
% imagesc
% jet
% caxis([0 max_value*0.7])
%
% 不进行Gamma增强
% 不进行分层归一化
% 不进行背景扣除
%% ================================================================

figure('Color','w', ...
    'Position',[150 80 900 750]);

imagesc(TFM_env);

set(gca,'YDir','normal');

axis equal tight;

colormap(jet);

caxis([0 max_value*color_ratio]);

colorbar;

xlabel('X (mm)');
ylabel('Depth (mm)');

title('Multi-layer TFM Reconstruction');

hold on;

if show_interface

    plot([1 Nx],[z_int1 z_int1], ...
        'w--','LineWidth',1.2);

    plot([1 Nx],[z_int2 z_int2], ...
        'w--','LineWidth',1.2);

end

hold off;

%% ================================================================
% 13. 保存结果
%% ================================================================

if save_result

    save_path = fullfile( ...
        file_path,'TFM_original_lowmemory.mat');

    save(save_path, ...
        'TFM_img', ...
        'TFM_env', ...
        't0_index', ...
        'probeX', ...
        'dt', ...
        'h1','h2','h3', ...
        'v1','v2','v3', ...
        '-v7');

    fprintf('\n数据保存完成：\n%s\n',save_path);

end

fprintf('\n========================================\n');
fprintf('          TFM反演全部完成\n');
fprintf('========================================\n');

%% ================================================================
% 局部函数：原始线性插值
%
% 等价于原代码：
%
% it0 = floor(it);
% alpha = it-it0;
%
% S(valid) =
% (1-alpha).*sig(idx)+alpha.*sig(idx+1);
%% ================================================================

function S = tfmInterp(sig,it)

    Nt = numel(sig);

    S = zeros(size(it),'single');

    it0 = floor(it);

    alpha = it-it0;

    valid = ...
        (it0>=1) & ...
        (it0<Nt);

    if any(valid(:))

        idx = double(it0(valid));

        S(valid) = ...
            (1-alpha(valid)).*sig(idx) + ...
            alpha(valid).*sig(idx+1);

    end

end
