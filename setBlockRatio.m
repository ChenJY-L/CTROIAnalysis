function [subMasks, ratio_u, ratio_v] = setBlockRatio(firstImage, roiData)

% --- 1. 获取并校验用户输入的分割比例 ---
inputRatioPrompt = {'请输入沿P1-P2和P4-P3边的分割比例 (例如: [1 2 1])', ...
    '请输入沿P1-P4和P2-P3边的分割比例 (例如: [1 1])'};
defaultRatios = {'[1 1]', '[1 1]'}; % 默认值
inputTitle = '设置ROI分割比例';

inputResult = inputdlg(inputRatioPrompt, inputTitle, [1 50], defaultRatios);

if isempty(inputResult)
    disp('用户取消了输入，ROI分割未执行。');
    return; % 或者根据需要进行其他处理
end

try
    ratio_u = str2num(inputResult{1}); % u方向的比例 (对应 P1-P2, P4-P3)
    ratio_v = str2num(inputResult{2}); % v方向的比例 (对应 P1-P4, P2-P3)

    if isempty(ratio_u) || isempty(ratio_v) || ...
            ~isvector(ratio_u) || ~isvector(ratio_v) || ...
            any(ratio_u <= 0) || any(ratio_v <= 0) || ...
            ~isnumeric(ratio_u) || ~isnumeric(ratio_v)
        errordlg('分割比例输入无效。请输入由正数组成的向量，例如 "[1 2 1]"。', '输入错误');
        return;
    end
catch ME
    errordlg(['解析分割比例时出错: ', ME.message], '输入错误');
    return;
end

% --- 2. 提取ROI的四个角点 (假设顺序为P1, P2, P3, P4) ---
% P1 --u--> P2
%  ^        ^
%  v        v
% P4 --u--> P3
% (P1对应uv参数空间的(0,0), P2(1,0), P3(1,1), P4(0,1))
if size(roiData, 1) ~= 4
    errordlg('roiData必须包含4个坐标点才能进行四边形分割。', '数据错误');
    return;
end
P1 = roiData(1, :);
P2 = roiData(2, :);
P3 = roiData(3, :); % 注意：P3是(u=1,v=1)的点
P4 = roiData(4, :); % 注意：P4是(u=0,v=1)的点

% --- 3. 定义双线性插值函数 ---
% P(u,v) = (1-u)(1-v)P1 + u(1-v)P2 + u*v*P3 + (1-u)v*P4
bilinear_interp = @(u, v) (1-u).*(1-v).*P1 + u.*(1-v).*P2 + u.*v.*P3 + (1-u).*v.*P4;

% --- 4. 计算归一化的分割点坐标 (u和v方向) ---
% u_coords 将是 [0, r1_1/sum(r1), (r1_1+r1_2)/sum(r1), ..., 1]
u_coords = [0, cumsum(ratio_u) / sum(ratio_u)];
v_coords = [0, cumsum(ratio_v) / sum(ratio_v)];

num_u_segments = length(ratio_u); % u方向的子ROI数量
num_v_segments = length(ratio_v); % v方向的子ROI数量

subMasks = cell(1, num_u_segments * num_v_segments);
mask_idx = 0; % 使用从0开始的索引，方便计算

% --- 5. 循环生成每个子ROI的顶点 ---
for j = 1:num_v_segments      % 遍历v方向的每一段 (行)
    v_start = v_coords(j);    % 当前子ROI的起始v参数
    v_end   = v_coords(j+1);  % 当前子ROI的结束v参数

    for i = 1:num_u_segments  % 遍历u方向的每一段 (列)
        u_start = u_coords(i);    % 当前子ROI的起始u参数
        u_end   = u_coords(i+1);  % 当前子ROI的结束u参数

        mask_idx = mask_idx + 1;

        % 计算子ROI的四个角点 (顺序: 左下, 右下, 右上, 左上 - 对应参数空间)
        % P_sub1: (u_start, v_start)
        % P_sub2: (u_end,   v_start)
        % P_sub3: (u_end,   v_end)
        % P_sub4: (u_start, v_end)
        p_sub1 = bilinear_interp(u_start, v_start);
        p_sub2 = bilinear_interp(u_end,   v_start);
        p_sub3 = bilinear_interp(u_end,   v_end);
        p_sub4 = bilinear_interp(u_start, v_end);

        subMasks{mask_idx} = [p_sub1; p_sub2; p_sub3; p_sub4];
    end
end

disp(['成功生成 ', num2str(length(subMasks)), ' 个子ROI。']);

% --- 可选: 显示子ROI以供验证 ---
figure;
imshow(firstImage)
hold on;
% 绘制原始ROI
patch(roiData(:,1), roiData(:,2), 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Original ROI');

colors = jet(length(subMasks)); % 为每个子ROI使用不同颜色
for k = 1:length(subMasks)
    sub_roi_points = subMasks{k};
    if ~isempty(sub_roi_points) && size(sub_roi_points,2) >= 2 % 确保是2D或3D点
        patch(sub_roi_points(:,1), sub_roi_points(:,2), colors(k,:), ...
            'FaceAlpha', 0.5, 'EdgeColor', 'k', ...
            'DisplayName', ['Sub-ROI ' num2str(k)]);
        % 计算子ROI中心并标注数字 (仅2D)
        if size(sub_roi_points,2) == 2
            center_x = mean(sub_roi_points(:,1));
            center_y = mean(sub_roi_points(:,2));
            text(center_x, center_y, num2str(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end
axis equal;
title('分割后的子ROI');
legend show;
hold off;

end