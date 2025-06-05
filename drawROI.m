function [roiMask, roiData, roiShapeChoice] = drawROI(dicomDir, firstImage, unitX, unitY)
% Display the first image and let user draw ROI
figure('Name', '在图像上绘制ROI');
imshow(firstImage, []); % Display using full data range
title('请在图像上用鼠标绘制感兴趣区域 (ROI)。');

roiShapePrompt = {'请选择ROI的形状:\n1: 矩形 (Rectangle)\n2: 椭圆 (Ellipse)\n3: 多边形 (Polygon)\n4: 自由手 (Freehand)\n5: 从文件加载'};
roiShapeDlgtitle = '选择ROI形状';
roiShapeDefinput = {'3'}; % Default freehand
roiShapeResponse = inputdlg(roiShapePrompt, roiShapeDlgtitle, [1 40], roiShapeDefinput);
if isempty(roiShapeResponse)
    disp('ROI形状选择取消。程序退出。');
    close(gcf); % Close image window
    return;
end

roiShapeChoice = str2double(roiShapeResponse{1});
hROI = []; % Handle for the ROI graphics object

switch roiShapeChoice
    case 1 % Rectangle
        disp('请在图像上绘制一个矩形ROI (拖动绘制)。');
        hROI = drawrectangle;
    case 2 % Ellipse
         disp('请在图像上绘制一个椭圆ROI (拖动绘制)。');
         hROI = drawellipse;
    case 3 % Polygon with distance display
         disp('请在图像上绘制一个多边形ROI (双击结束绘制)。');
         hROI = drawPolygonWithDistance(unitX, unitY);
    case 4 % Freehand
         disp('请在图像上绘制一个自由手ROI (双击结束绘制)。');
         hROI = drawfreehand;
    case 5 % Load
         disp('请选择ROI的MAT文件')
         [matFileName, matFilePath] = uigetfile('*.mat', '选择ROI掩码文件 (例如: roiMask.mat)');
         
    otherwise
        disp('无效的ROI形状选择。程序退出。');
        close(gcf); % Close image window
        return;
end

% Wait for user to finish drawing the ROI
if ~isempty(hROI) & roiShapeChoice~=5
    if roiShapeChoice ~= 3  % For non-polygon shapes, use wait
        wait(hROI);
    end
end

% Check if a valid ROI object was drawn
if (isempty(hROI) || ~isvalid(hROI)) && roiShapeChoice~=5
     disp('ROI绘制失败或取消。程序退出。');
     close(gcf); % Close image window
     return;
end

% Create the binary mask from the drawn ROI
if roiShapeChoice ~= 5
    roiMask = createMask(hROI); % roiMask will be MxN
    roiData = hROI.Position;
else
    loadedData = load(fullfile(matFilePath, matFileName));
    roiMask = loadedData.roiMask;
    roiData = loadedData.roiData;
    roiShapeChoice = loadedData.roiShapeChoice;
end

% Close the ROI drawing figure
saveas(gcf, fullfile(dicomDir, "ROI.png"))
save(fullfile(dicomDir, "ROI.mat"), "roiMask", "roiData", "roiShapeChoice")
close(gcf);

% Check if the drawn ROI area is empty
if ~any(roiMask(:))
    disp('绘制的ROI区域为空。程序退出。');
    return;
end

%
fprintf('ROI绘制完成。\n');
end