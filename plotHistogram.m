function plotHistogram(allHistogramsCell, allTimeStampCell, histBinEdges, threshold, eventDateTimes, eventLabels)

histBinCenters = (histBinEdges(1:end - 1) + histBinEdges(2:end)) / 2;

allHistograms = [];
allTimeStamp = [];
for i = 1:numel(allHistogramsCell)
    allHistograms = [allHistograms; allHistogramsCell{i}];
    allTimeStamp = [allTimeStamp, allTimeStampCell{i}];
end

figure('Name', '直方图时间演变热图', 'Position', [100, 100, 1400, 800]);
tiledlayout('flow')
nexttile,
% 创建时间索引
timeIndices = 1:size(allHistograms, 1);

%% 绘制热图
imagesc(histBinCenters, timeIndices, allHistograms);
colorbar;
xlabel('像素值');
ylabel('帧索引');
title('ROI内像素值分布随时间的变化 (热图)');
colormap('hot');

% 添加阈值线
hold on;
plot([threshold, threshold], [1, length(timeIndices)], 'g--', 'LineWidth', 2, 'DisplayName', '阈值');
legend('Location', 'best');
hold off;
%%
nexttile,
threshHotMap = cumsum(allHistograms, 2);

% 绘制热图
imagesc(histBinCenters, timeIndices, threshHotMap);
colorbar;
xlabel('像素值');
ylabel('帧索引');
title('ROI内像素值分布随时间的变化 (统计图)');
colormap('hot');

%%
nexttile,
% 选择几个代表性时间点
diffTimeStamp = abs(allTimeStamp(:)' - eventDateTimes(:));
[~, selectedFrames] = min(diffTimeStamp, [], 2);
colors = lines(length(selectedFrames));
hold on;

for i = 1:length(selectedFrames)
    frameIdx = selectedFrames(i);
    plot(histBinCenters, allHistograms(frameIdx, :), ...
         'Color', colors(i,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('%s', eventLabels{i}));
end

% 添加阈值线
ylim_current = ylim;
plot([threshold, threshold], ylim_current, 'k--', 'LineWidth', 2, 'DisplayName', '阈值');

xlabel('像素值');
ylabel('归一化频率');
title('ROI内不同时间点的像素值分布对比');
legend('Location', 'best');
grid on;
hold off;

%%
nexttile,
hold on
for i = 1:length(selectedFrames)
    frameIdx = selectedFrames(i);
    plot(histBinCenters, threshHotMap(frameIdx, :), ...
         'Color', colors(i,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('%s', eventLabels{i}));
end

xlabel('像素值');
ylabel('归一化频率');
title('ROI内不同时间点的像素值分布对比');
legend('Location', 'best');
grid on;
hold off;

end