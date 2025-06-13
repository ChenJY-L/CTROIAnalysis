function plotSubROI(allTimeStampCell, allCountsBelowCellSub, subMasks, ...
                    ratio_u, ratio_v, ...
                    eventTimesStr, eventLabels, threshold)
figure('Name', 'SubROI内像素计数变化曲线 (按帧)', 'Position', [500, 500, 1600, 500]);

tiledlayout(length(ratio_u), length(ratio_v)), 
% --- Define colors and markers ---
num_u = length(ratio_u);
plotColors = lines(max(1, length(ratio_v))); % One color per row
plotMarkers = {'o', 's', '^', 'd', 'v', 'p', 'h', 'x', '*', '+'}; % Markers for columns
lineStyles = {'-', '--', ':', '-.'}; % Line styles can also be varied
colors = parula(num_u);

% 将字符串时间转换为 datetime 对象
eventDateTimes = datetime(eventTimesStr, 'InputFormat', 'yyyyMMdd HH:mm');

for s_idx = 1:length(subMasks) 
    nexttile, 
    hold on
    for i = 1:numel(allTimeStampCell)
        markIndex = mod(s_idx - 1, num_u) + 1;
        % colorIndex = floor((s_idx - 1) / num_u) + 1;
        
        plot(allTimeStampCell{i}, allCountsBelowCellSub{i}(s_idx, :)*100, ...
            'LineStyle', lineStyles{markIndex}, ...
            'Marker', plotMarkers{markIndex}, ...
            'Color', [0, 0.4470, 0.7410]...
            );
    end

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
    title(sprintf('ROI区域内低于阈值%g的像素数量占比 -- Sub%d', threshold, s_idx));
    % legend('Location', 'best');
    grid on;
end


end