function plotResult(allTimeStampCell, ...
    allCountsBelowCell, allEdgePositionCell, allMeanCell, allHistogramsCell, ...
    histBinEdges, ...
    eventTimesStr, eventLabels, threshold, unitY)

figure('Name', 'ROI内像素计数变化曲线 (按帧)', 'Position', [500, 500, 1600, 500]);
set(0, 'DefaultAxesFontSize', 14); % 设置坐标轴字体大小
set(0, 'DefaultTextFontSize', 14); % 设置文本字体大小
tiledlayout('flow')
% 将字符串时间转换为 datetime 对象
eventDateTimes = datetime(eventTimesStr, 'InputFormat', 'yyyyMMdd HH:mm');

%% Plot counts against cumulative frame index
nexttile, 
hold on; % Keep current plot
meanTime = {};
meanCounts = [];
for i = 1:numel(allTimeStampCell)
    meanTime{i} = allTimeStampCell{i}(1);
    meanCounts(i) = mean(allCountsBelowCell{i}*100);
    isOut = isoutlier(allCountsBelowCell{i}, 'quartiles'); 

    plot(allTimeStampCell{i}(~isOut), allCountsBelowCell{i}(~isOut)*100, ...
        'o',  ...
        'Color', [0, 0.4470, 0.7410]);
end
plot(datetime([meanTime{:}]), meanCounts, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]) 
% 添加 xlines
for k = 1:numel(eventDateTimes)
    h = xline(eventDateTimes(k), 'r--', eventLabels{k});
    h.FontSize = 13;
end
hold off; % Release plot

% Add labels, title, legend
xlabel("时间");
ylabel('像素数量(归一化)');
ytickformat('percentage')
title(sprintf('ROI区域内低于阈值%g的像素数量占比', threshold), 'FontSize', 18);
% legend('Location', 'best');
grid on;
box on;

%% 
nexttile, 
meanTime = {};
meanEdge = [];
hold on;
for i = 1:numel(allTimeStampCell)
    meanTime{i} = allTimeStampCell{i}(1);
    meanEdge(i) = mean(allEdgePositionCell{i})*unitY;
    isOut = isoutlier(allEdgePositionCell{i}, 'quartiles'); 
    plot(allTimeStampCell{i}(~isOut), allEdgePositionCell{i}(~isOut).*unitY, 's', 'Color', [0.8500, 0.3250, 0.0980]);
end
set(gca, 'YDir', 'reverse');
plot(datetime([meanTime{:}]), meanEdge, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]) 

% 添加事件标记
for k = 1:numel(eventDateTimes)
    h = xline(eventDateTimes(k), 'r--', eventLabels{k});
    h.FontSize = 13;
end
hold off;

xlabel("时间");
ylabel('表皮层平均像素y坐标(mm)');
title('表皮层位移', 'FontSize', 18);
grid on;
box on;
%% 
nexttile
 
meanTime = {};
meanMean = [];
hold on;
for i = 1:numel(allTimeStampCell)
    meanTime{i} = allTimeStampCell{i}(1);
    meanMean(i) = mean(allMeanCell{i});
    isOut = isoutlier(allMeanCell{i}, 'quartiles'); 
    plot(allTimeStampCell{i}(~isOut), allMeanCell{i}(~isOut), ...
        '-o',  ...
        'Color', [0.4940 0.1840 0.5560]);
end
plot(datetime([meanTime{:}]), meanMean, '--', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840]) 

% 添加 xlines
for k = 1:numel(eventDateTimes)
    h = xline(eventDateTimes(k), 'r--', eventLabels{k});
    h.FontSize = 13;
end
hold off; % Release plot

% Add labels, title, legend
xlabel("时间");
ylabel('灰度均值');
grid on;
box on;
title('灰度均值')

%% 设置base
baseIdx = inputdlg({"输入Base的索引"}, ...
                    "Base的索引", [1 35], {"1"});
baseIdx = str2double(baseIdx{1});
%% 绘制直方图
nexttile,
histBinCenters = (histBinEdges(1:end - 1) + histBinEdges(2:end)) / 2;
% hold on;
% for i = 1:numel(allTimeStampCell)
%     [~, idx] = min(abs(allCountsBelowCell{i} - meanCounts(i)));
% 
%     plot(histBinCenters, allHistogramsCell{i}(idx, :), ...
%         '-', 'DisplayName', datestr(allTimeStampCell{i}(idx), 'HH:MM:SS')); 
% end
% hold off;
% legend
% xlabel('灰度值')
% ylabel('像素数量')


frameTime = zeros(1, numel(allTimeStampCell));
histogramData = zeros(numel(allTimeStampCell), length(histBinCenters));
for i = 1:numel(allTimeStampCell)
    [~, idx] = min(abs(allCountsBelowCell{i} - meanCounts(i)));

    frameTime(i) = datenum(allTimeStampCell{i}(idx));
    histogramData(i, :) = allHistogramsCell{i}(idx, :);
end
histogramData = histogramData - histogramData(baseIdx, :);

% waterfall(histBinCenters, frameTime(baseIdx:end), histogramData(baseIdx:end, :));
% plot(histBinCenters, histogramData(baseIdx:end, :))
% legend(datestr(frameTime(baseIdx:end), 'HH:MM'))

imagesc(histBinCenters, frameTime, histogramData)
hold on
for k = 1:numel(eventDateTimes)
    plot([0, 255], ...
         [datenum(eventDateTimes(k)), datenum(eventDateTimes(k))], ...
          '--', 'DisplayName', eventLabels{k});
end
legend
hold off
colormap('jet')
xlabel('灰度值')
ylabel('时间')
datetick('y', 'HH:MM:SS', 'keepticks'); 
title('直方图变化量')
colorbar

%% 绘制各个分量的变化趋势
nexttile,
% dataMatrix = zeros(length() ,length(histBinCenters))
dataMatrix = allHistogramsCell{1};
dataTime = allTimeStampCell{1};
dataEdge = allEdgePositionCell{1};
dataMean = allMeanCell{1};
dataMeanMatrix = mean(allHistogramsCell{1});
for i = 2:numel(allTimeStampCell) 
    dataTime = [dataTime, allTimeStampCell{i}];
    dataMatrix = [dataMatrix; allHistogramsCell{i}];
    dataEdge = [dataEdge, allEdgePositionCell{i}];
    dataMean = [dataMean, allMeanCell{i}];

    isOut = isoutlier(allEdgePositionCell{i}, 'quartiles'); 
    dataMeanMatrix = [dataMeanMatrix; mean(allHistogramsCell{i}(~isOut, :))];
end

baseData = mean(allHistogramsCell{baseIdx}, 1);
diffMeanMatrix = dataMeanMatrix - baseData;
diffMatrix = dataMatrix - baseData;

hold on 
% plot(dataTime, diffMatrix)
plot(datetime([meanTime{:}]), diffMeanMatrix, '-.');
hold off

xlabel("时间");
ylabel("像素数量（归一化）")
grid on;
box on;

% --- 添加图例 ---
numBins = size(dataMatrix, 2);
% legendLabels = cell(1, numBins);
% for k = 1:numBins
%     legendLabels{k} = sprintf('灰度值 %.0f ~ %.0f', histBinEdges(k), histBinEdges(k + 1));
% end
% legend(legendLabels, 'Location', 'eastoutside', 'FontSize', 8, 'NumColumns', 2); % 将图例放在图的外部，防止遮挡

%% 是否保存数据
saveChoice = questdlg('是否要保存分析数据 (原始数据和差值数据)?', ...
    '保存数据', ...
    '是', '否', '否');

if strcmp(saveChoice, '是')
    % 弹出文件保存对话框
    [file, path] = uiputfile('*.xlsx', '选择保存位置和文件名', 'ProcessedData.xlsx');
    
    % 如果用户确认保存
    if ischar(file)
        fullFileName = fullfile(path, file);
        
        % 1. 创建所有表格都共用的列 (只执行一次)
        commonVars = table(dataTime', dataEdge'.*unitY, dataMean', ...
            'VariableNames', {'时间', '表皮层位移_mm', '灰度均值'});
        
        % 2. 创建灰度区间的列名模板 (只执行一次)
        numBins = size(dataMatrix, 2);
        binLabels = arrayfun(@(k) sprintf('Gray_%.0f_%.0f', histBinEdges(k), histBinEdges(k+1)), ...
                             1:numBins, 'UniformOutput', false);

        % 3. 创建并写入第一个Sheet：'差值数据'
        T_diff_main = array2table(diffMatrix, 'VariableNames', strcat('Diff_', binLabels));
        writetable([commonVars, T_diff_main], fullFileName, 'Sheet', '差值数据');

        % 4. 创建并写入第二个Sheet：'原始数据'
        T_raw_main = array2table(dataMatrix, 'VariableNames', strcat('Raw_', binLabels));
        writetable([commonVars, T_raw_main], fullFileName, 'Sheet', '原始数据');
        
        % 5. 创建第三个Sheet
        T_mean_main = array2table(diffMeanMatrix, 'VariableNames', strcat('Mean_Diff_', binLabels));
        T_mean_main = addvars(T_mean_main, datetime([meanTime{:}])', 'Before', 1, 'NewVariableNames', '时间');
        writetable(T_mean_main, fullFileName, 'Sheet', '均值变化量');

        T_mean_main = array2table(dataMeanMatrix, 'VariableNames', strcat('Mean_', binLabels));
        T_mean_main = addvars(T_mean_main, datetime([meanTime{:}])', 'Before', 1, 'NewVariableNames', '时间');
        writetable(T_mean_main, fullFileName, 'Sheet', '均值数据');
        
        % 6. 完成提示
        % msgbox(sprintf('数据已成功保存到文件:\n%s', fullFileName), '保存成功');
        fprintf('所有数据已成功保存到: %s\n', fullFileName);
        
    else
        disp('用户取消了保存操作。');
    end
end

end