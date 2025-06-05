%% Custom function for drawing polygon with distance display
function hROI = drawPolygonWithDistance(unitX, unitY)
% Get current axes
ax = gca;

% Initialize variables
vertices = [];
lastPoint = [];
lastlastPoint = [];
lineHandles = [];
distanceText = [];
tempLineHandle = [];
isPolygonClosed = false; % Flag to track if polygon is closed

% Set up mouse callbacks
set(gcf, 'WindowButtonDownFcn', @mouseClick);
set(gcf, 'WindowButtonMotionFcn', @mouseMove);

% Instructions
title('单击添加多边形顶点，双击完成绘制。移动鼠标查看到上一个点的距离(mm)。');

% Wait for polygon completion
uiwait(gcf);

% Create the polygon ROI object
if size(vertices, 1) >= 3
    hROI = drawpolygon('Position', vertices);
else
    hROI = [];
    warning('多边形至少需要3个点。');
end

% Clean up
if ~isempty(distanceText) && isvalid(distanceText)
    delete(distanceText);
end
if ~isempty(tempLineHandle) && isvalid(tempLineHandle)
    delete(tempLineHandle);
end

function mouseClick(~, ~)
    % Get current point
    cp = get(ax, 'CurrentPoint');
    x = cp(1,1);
    y = cp(1,2);

    % Check if click is within axes limits
    xlim = get(ax, 'XLim');
    ylim = get(ax, 'YLim');
    if x < xlim(1) || x > xlim(2) || y < ylim(1) || y > ylim(2)
        return;
    end

    % Check for double click
    if strcmp(get(gcf, 'SelectionType'), 'alt')
        % Double click - finish polygon
        if size(vertices, 1) >= 3
            % Close the polygon by connecting to first point
            hold on;
            plot(ax, [vertices(end,1), vertices(1,1)], [vertices(end,2), vertices(1,2)], 'r-', 'LineWidth', 2);
            hold off;

            % Set polygon as closed and clean up distance display
            isPolygonClosed = true;
            if ~isempty(distanceText) && isvalid(distanceText)
                delete(distanceText);
            end
            if ~isempty(tempLineHandle) && isvalid(tempLineHandle)
                delete(tempLineHandle);
            end

            % Update title to indicate completion
            title('多边形ROI绘制完成');

            uiresume(gcf);
        else
            disp('多边形至少需要3个点。请继续添加点。');
        end
        return;
    end

    if ~isPolygonClosed
        % Single click - add vertex
        vertices = [vertices; x, y];

        % Plot the point
        hold on;
        plot(ax, x, y, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

        % Update last point to current point
        if ~isempty(lastPoint)
            lastlastPoint = lastPoint;
        end
        lastPoint = [x, y];

        % Draw line to previous point (if exists)
        if size(vertices, 1) > 1
            prevVertex = vertices(end-1, :);
            lineHandle = plot(ax, [prevVertex(1), x], [prevVertex(2), y], 'r-', 'LineWidth', 2);
            lineHandles = [lineHandles, lineHandle];
        end

        hold off;
    end
end

function mouseMove(~, ~)
    % Don't show distance if polygon is closed
    if isPolygonClosed
        return;
    end

    % Only show distance after first point is set
    if isempty(lastPoint)
        return;
    end

    % Get current mouse position
    cp = get(ax, 'CurrentPoint');
    x = cp(1,1);
    y = cp(1,2);

    % Check if mouse is within axes limits
    xlim = get(ax, 'XLim');
    ylim = get(ax, 'YLim');
    if x < xlim(1) || x > xlim(2) || y < ylim(1) || y > ylim(2)
        return;
    end

    % Calculate distance to last point in pixels
    distancePixels = sqrt((x - lastPoint(1))^2 + (y - lastPoint(2))^2);

    % Convert to mm using pixel spacing
    % Use average of unitX and unitY for distance calculation
    avgPixelSpacing = (unitX + unitY) / 2;
    distanceMM = distancePixels * avgPixelSpacing;
    
    theta = 0;
    if ~isempty(lastlastPoint)
        A = [x,y] - lastPoint;
        B = lastPoint - lastlastPoint;
        theta = rad2deg(acos(dot(A,B)/(norm(A)*norm(B))));
    end

    % Delete previous distance display
    if ~isempty(distanceText) && isvalid(distanceText)
        delete(distanceText);
    end
    if ~isempty(tempLineHandle) && isvalid(tempLineHandle)
        delete(tempLineHandle);
    end

    % Draw temporary line to last point
    hold on;
    tempLineHandle = plot(ax, [x, lastPoint(1)], [y, lastPoint(2)], 'g--', 'LineWidth', 1);
    
    % Display distance text in mm
    midX = (x + lastPoint(1)) / 2;
    midY = (y + lastPoint(2)) / 2;
    distanceText = text(midX, midY, sprintf('%.2f mm, %.2f° ', distanceMM, theta), ...
        'Color', 'green', 'FontSize', 10, 'FontWeight', 'bold', ...
        'BackgroundColor', 'white', 'EdgeColor', 'green', ...
        'HorizontalAlignment', 'center');

    hold off;
end
end