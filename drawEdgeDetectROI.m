function edgeRoiMask = drawEdgeDetectROI(firstImage, dicomDir)
% 5.5 绘制边缘检测专用ROI
figure('Name', '在图像上绘制边缘检测ROI');
imshow(firstImage, []); 
title('请在图像上绘制用于边缘检测的感兴趣区域 (Edge ROI)。');

% 获取边缘检测ROI形状选择
edgeRoiShapePrompt = {'请选择边缘检测ROI的形状:\n1: 矩形 (Rectangle)\n2: 椭圆 (Ellipse)\n3: 多边形 (Polygon)\n4: 自由手 (Freehand)\n5:从文件加载'};
edgeRoiShapeDlgtitle = '选择边缘检测ROI形状';
edgeRoiShapeDefinput = {'1'}; % 默认矩形
edgeRoiShapeResponse = inputdlg(edgeRoiShapePrompt, edgeRoiShapeDlgtitle, [1 40], edgeRoiShapeDefinput);

if isempty(edgeRoiShapeResponse)
    disp('边缘检测ROI形状选择取消。程序退出。');
    close(gcf);
    return;
end

edgeRoiShapeChoice = str2double(edgeRoiShapeResponse{1});
hEdgeROI = [];

switch edgeRoiShapeChoice
    case 1 % Rectangle
        disp('请在图像上绘制一个矩形边缘检测ROI。');
        hEdgeROI = drawrectangle;
    case 2 % Ellipse
        disp('请在图像上绘制一个椭圆边缘检测ROI。');
        hEdgeROI = drawellipse;
    case 3 % Polygon
        disp('请在图像上绘制一个多边形边缘检测ROI。');
        hEdgeROI = drawpolygon;
    case 4 % Freehand
        disp('请在图像上绘制一个自由手边缘检测ROI。');
        hEdgeROI = drawfreehand;
    case 5
        disp('请选择ROI的MAT文件')
        [matFileName, matFilePath] = uigetfile('*.mat', '选择ROI掩码文件 (例如: roiMask.mat)');
    otherwise
        disp('无效的边缘检测ROI形状选择。程序退出。');
        close(gcf);
        return;
end

% 创建边缘检测ROI掩码
if edgeRoiShapeChoice~=5
    % 等待用户完成绘制
    if ~isempty(hEdgeROI)
        wait(hEdgeROI);
    end

    if ~isempty(hEdgeROI) && isvalid(hEdgeROI) 
        edgeRoiMask = createMask(hEdgeROI);
    else
        disp('边缘检测ROI绘制失败。程序退出。');
        close(gcf);
        return;
    end
else
    loadedData = load(fullfile(matFilePath, matFileName));
    edgeRoiMask = loadedData.edgeRoiMask;
end

% 保存边缘检测ROI
saveas(gcf, fullfile(dicomDir, "EdgeROI.png"));
save(fullfile(dicomDir, "EdgeROI.mat"), "edgeRoiMask");
close(gcf);

fprintf('边缘检测ROI绘制完成。\n');
end