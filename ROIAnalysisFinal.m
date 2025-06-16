% Clear workspace, close figures, clear command window 
clear;
clc;
% close all;

%% 1. Select the directory containing DICOM files
dicomDir = uigetdir('', '请选择包含DICOM文件的文件夹');
if dicomDir == 0
    disp('用户取消了文件夹选择。程序退出。');
    return; % Exit the script
end

%% 2. List all files in the directory and identify potential DICOMs
[dicomFiles, filePaths] = scanDicomFiles(dicomDir);

numFiles = length(dicomFiles);
if numFiles == 0
    disp('在指定的文件夹中没有找到有效的DICOM文件。程序退出。');
    return;
else
    fprintf('找到 %d 个有效的DICOM文件。\n', numFiles);
end

%% 3. Get user input for number of files to skip
filesToProcessIdx = setSkipFiles(numFiles);

numFilesToProcess = length(filesToProcessIdx);
if numFilesToProcess == 0
     disp('根据跳过设置，没有文件需要处理。程序退出。');
     return;
end

%% 4. Read the first FRAME of the FIRST DICOM file (original index) to define ROI
firstFilePath = filePaths{filesToProcessIdx(1)};
try
    % Read DICOM info to get pixel spacing
    imageInfo = dicominfo(firstFilePath);
    unitX = imageInfo.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX*10; % mm per pixel in X direction
    unitY = imageInfo.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY*10; % mm per pixel in Y direction
    
    % Read the entire volume/multiframe data from the first file
    firstFileVolume = dicomread(firstFilePath, "Frames", 1);
    % Extract the first frame (assuming frame dimension is 3rd or higher)
    volSize = size(firstFileVolume);

    firstImage = firstFileVolume(:,:,1); % Assume 3rd dim is frame/slice
catch ME
    error(['读取第一个DICOM文件或其第一帧失败，无法绘制ROI: ', firstFilePath, newline, ME.message]);
end


%% 5. Get user selected ROI shape and draw
[roiMask, roiData, roiShapeChoice] = drawROI(dicomDir, firstImage, unitX, unitY);

%% 5.5 将ROI分块处理
useSubROIs = false;
if roiShapeChoice == 3 
    useSubROIs = true;
    [subMasks, ratio_u, ratio_v] = setBlockRatio(firstImage, roiData);
end

% 设置边缘检测ROI
edgeRoiMask = drawEdgeDetectROI(firstImage, dicomDir);

%% 6. Get user input for the threshold value
promptThreshold = {'请输入像素值阈值:'};
dlgtitleThreshold = '设置阈值';
definputThreshold = {'30'}; % Default threshold example
inputResponseThreshold = inputdlg(promptThreshold, dlgtitleThreshold, [1 35], definputThreshold);

if isempty(inputResponseThreshold)
    disp('阈值输入取消。程序退出。');
    return;
end

threshold = str2double(inputResponseThreshold{1});

if isnan(threshold)
    disp('输入的阈值无效。程序退出。');
    return;
end

fprintf('阈值设置为: %g\n', threshold);

%% 7. Initialize arrays to store results (dynamically grow)
allCountsBelowCell = cell(1, length(filesToProcessIdx));
allTimeStampCell = cell(1, length(filesToProcessIdx));
% 添加边缘检测结果存储
allEdgePositionCell = cell(1, length(filesToProcessIdx));
allCountsBelowCellSub = cell(1, length(filesToProcessIdx));
allHistogramsCell = cell(1, length(filesToProcessIdx)); % 存储每个文件的直方图数据
allElipseCoeffs = cell(1, length(filesToProcessIdx));
allMeanCell = cell(1, length(filesToProcessIdx));
processedFrameCount = 0;
refEdgePosition = zeros(0);

%% 8. Loop through the selected DICOM files, process frames, apply ROI/threshold, and count
fprintf('正在处理文件和其中的帧...\n');
ROIArea = sum(roiMask, 'all');
frameInterval = inputdlg("帧间隔", "帧间隔设置", [1 35], "15");
frameInterval = str2double(frameInterval{1});

histBinEdges = linspace(0, 255, 51 + 1);
histBinEdges(2:end) = histBinEdges(2:end) + 0.1;

corrFlag = inputdlg("是否根据表皮层位移对ROI进行矫正", "ROI矫正设置", [1 35], "1");
corrFlag = logical(str2double(corrFlag));

saveVideoFlag = inputdlg("是否保存为视频", "Save", [1 35], "0");
saveVideoFlag = logical(str2double(saveVideoFlag));

grayCorr = inputdlg({"是否进行灰度值矫正", "矫正斜率K(1/mm)"}, ...
                        "灰度矫正", [1 35], {"0", "13.513"});
grayCorrFlag = logical(str2double(grayCorr{1}));

if grayCorrFlag
    grayCorrK = str2double(grayCorr{2});
    % grayCorrb = str2double(grayCorr{3});
    grayCorrFunc = @(h) grayCorrK*h;
end

%% 8.5 循环处理
if saveVideoFlag
    figure('Name', 'CurrentFrame')
    videoFile = fullfile(dicomDir,'output.avi');  % 或 'output.mp4'
    v = VideoWriter(videoFile);  % 默认 AVI 格式
    v.FrameRate = 5;             % 可自定义帧率
    open(v);
end

% 保留初始ROI的位置
oriROIMask = roiMask;
if useSubROIs
    oriSubMasks = subMasks;
end

for k = filesToProcessIdx
    % #--- 文件读取 ---#
    currentFileName = dicomFiles{k};
    currentFilePath = filePaths{k};

    % Read the entire volume/multiframe data from the current file
    dFile = dicomFile(currentFilePath);
    imageDataVolume = getPixelData(dFile);
    % imageDataVolume = dicomread(currentFilePath);
    imageInfo = dicominfo(currentFilePath);
    timeStr = [imageInfo.AcquisitionDate, ' ', imageInfo.AcquisitionTime];
    startTime = datetime(timeStr, 'InputFormat', 'yyyyMMdd HHmmss.SSS');
    
    if size(imageDataVolume, 4) == 1
        frameTime = milliseconds(0.0333);
    else
        frameTime = milliseconds(imageInfo.FrameTime);
    end

    % Get dimensions of the volume
    volSize = size(imageDataVolume);
    M = volSize(1); % Image height
    N = volSize(2); % Image width
    P = 1; % Initialize frame/slice count within this file

    if length(volSize) >= 3
        P = 1; % 3rd dimension is often frames/slices
    end

    % Handle 4D data (e.g., MxNxPxT) - treat 4th dim as extra frames
    if length(volSize) == 4
        P = volSize(4); % Total frames = slices * time points
        M = volSize(1); N = volSize(2); % Image dimensions
        % Need to reshape or index carefully below
    elseif length(volSize) > 4
        warning(['文件 ', currentFileName, ' 维度大于4D，可能无法正确处理。跳过。']);
        continue; % Skip this file
    end

    if M ~= size(roiMask, 1) || N ~= size(roiMask, 2)
        warning(['图像尺寸不匹配ROI掩码，跳过文件: ', currentFileName]);
        continue; % Skip this file
    end

    if P == 0
        warning(['文件 ', currentFileName, ' 没有可处理的帧/切片。跳过。']);
        continue; % Skip this file
    end

    % Loop through each frame/slice within the current volume
    tmpPointsBelow = zeros(0);
    tmpEdgePosition = zeros(0);
    tmpTimeStamp = NaT(0);
    tmpHistograms = [];
    tmpElipseCoeffs = [];
    tmpMean = [];
    if useSubROIs
        currentFileCountsSubROI = zeros(length(subMasks), length(1:frameInterval:P));
    end
    
    % #--- 文件处理 ---#
    count = 0;
    for frameIdx = 1:frameInterval:P
        count = count + 1;
        processedFrameCount = processedFrameCount + 1; % Increment total frame counter

        % Extract the current frame
        if length(volSize) == 4
            % For 4D, calculate slice and time index
            sliceIdx = mod(frameIdx - 1, volSize(3)) + 1;
            timeIdx = floor((frameIdx - 1) / volSize(3)) + 1;
            currentFrame = squeeze(imageDataVolume(:,:,1, timeIdx)); % 灰度图像，仅读取第一个通道
        else % 2D, 3D (MxNxP) or other cases
            currentFrame = squeeze(imageDataVolume(:,:,frameIdx));
        end
        
        % 表皮层位置检测
        [rows, cols] = detectEpidermis(currentFrame, edgeRoiMask, [0.15, 0.22]);
        
        if processedFrameCount == 1
            % 保留第一帧中表皮层的位置
            refEdgePosition = rows;
        end

        % 根据表皮层相对于上一次测量的位移，修正ROI的位置
        if corrFlag & (~isempty(refEdgePosition))
            % 计算表皮层位移
            epidermisDiff = round(mean(refEdgePosition) - mean(rows)); % 参考表皮层与当前表皮层位置做差

            % 修正ROI位置
            roiMask = fixROI(oriROIMask, epidermisDiff);

            if useSubROIs
                for subMaskIdx = 1:length(subMasks)
                    subMasks{subMaskIdx}(:, 2) = oriSubMasks{subMaskIdx}(:, 2) - epidermisDiff; 
                end
            end
        end
        
        % # --- 绘制结果 --- #
        if saveVideoFlag
            % 调试测试
            imshow(labeloverlay(currentFrame, roiMask)),
            hold on 
            plot(cols, rows, 'ro', 'MarkerSize', 3)  
            if useSubROIs
                colors = jet(length(subMasks));
                for s = 1:length(subMasks)
                    sub_roi_points = subMasks{s};
                    if ~isempty(sub_roi_points) && size(sub_roi_points,2) >= 2 % 确保是2D或3D点
                        patch(sub_roi_points(:,1), sub_roi_points(:,2), colors(s,:), ...
                              'FaceAlpha', 0.5, 'EdgeColor', 'k', ...
                              'DisplayName', ['Sub-ROI ' num2str(s)]);
                        % 计算子ROI中心并标注数字 (仅2D)
                        if size(sub_roi_points,2) == 2
                            center_x = mean(sub_roi_points(:,1));
                            center_y = mean(sub_roi_points(:,2));
                            text(center_x, center_y, num2str(s), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                        end
                    end
                end
            end
            hold off
            drawnow
    
            % 写入处理结果
            frame = getframe(gca);
            writeVideo(v, frame);
        end
        
        % 根据耦合剂厚度进行灰度矫正
        if grayCorrFlag & (~isempty(refEdgePosition))
            agentHeight = (rows - refEdgePosition)*unitY;
            currentFrame(roiMask) = currentFrame(roiMask) - mean(grayCorrFunc(agentHeight));
        end

        % Apply the ROI mask: get pixel values within the ROI
        frameInsideROI = currentFrame(roiMask);
        
        % 估计衰减系数
        % [elipseCoeffs, elipseStat] = elipseValueEstimation(currentFrame, roiMask, rows, unitY);
        
        % Count pixels above and below/equal to the threshold within the ROI
        countsBelowFrame = sum(frameInsideROI <= threshold)/ROIArea;

        % 提取直方图数据
        [histCounts, ~] = histcounts(frameInsideROI, histBinEdges);
        histCounts = histCounts / ROIArea; % Normalize to probability
        
        % 提取子ROI的数据
        if useSubROIs
            for s_idx = 1:length(subMasks)
                subPoints = subMasks{s_idx};

                subPoints(:,1) = max(1, min(subPoints(:,1), N)); % X-coords
                subPoints(:,2) = max(1, min(subPoints(:,2), M)); % Y-coords

                subRoiGeneratedMask = poly2mask(subPoints(:,1), subPoints(:,2), M, N);
                subROIArea = sum(subRoiGeneratedMask, 'all');

                frameInsideSubROI = currentFrame(subRoiGeneratedMask);
                countsBelowFrameSub = sum(frameInsideSubROI <= threshold) / subROIArea;
                currentFileCountsSubROI(s_idx, count) = countsBelowFrameSub;
            end
        end

        % Store the counts for this frame (dynamically grow arrays)
        tmpTimeStamp(end+1) = startTime + frameIdx*frameTime;
        tmpPointsBelow(end+1) = countsBelowFrame;
        tmpEdgePosition(end+1) = mean(rows);
        tmpHistograms = [tmpHistograms; histCounts];
        tmpMean(end+1) = mean(frameInsideROI);

        % tmpElipseCoeffs = [tmpElipseCoeffs; elipseCoeffs'];
        
    end % End of frame loop
    
    idx = k - filesToProcessIdx(1) + 1;
    allTimeStampCell{idx} = tmpTimeStamp;
    allCountsBelowCell{idx} = tmpPointsBelow;
    allEdgePositionCell{idx} = tmpEdgePosition;
    allHistogramsCell{idx} = tmpHistograms;
    allMeanCell{idx} = tmpMean;

    % allElipseCoeffs{idx} = tmpElipseCoeffs;

    if useSubROIs
        allCountsBelowCellSub{idx} = currentFileCountsSubROI;
    end
    
    % 按照视频矫正
    % refEdgePosition = tmpEdgePosition;
    % oriROIMask = roiMask;
    % oriSubMasks = subMasks;
        
    % Optional: display progress per file
    fprintf('已处理文件 %d / %d (%s) - 包含 %d 帧。\n', ...
        find(filesToProcessIdx == k), numFilesToProcess, currentFileName, P);

end % End of file loop

if saveVideoFlag 
    close(v); 
end
%% 读取时间点的文件
% 读取CSV文件
[eventTimesStr, eventLabels] = readEvents(dicomDir);

%% 9. Plot the results
plotResult(allTimeStampCell, ...
    allCountsBelowCell, ...
    allEdgePositionCell, ...
    allMeanCell, ...
    allHistogramsCell, ...
    histBinEdges, ...
    eventTimesStr, eventLabels, ...
    threshold, unitY)

%% 10. 绘制SubMask的结果
if useSubROIs
    plotSubROI(allTimeStampCell, allCountsBelowCellSub, subMasks,...
               ratio_u, ratio_v, ...
               eventTimesStr, eventLabels, threshold)
end
%% 11. 直方图可视化
plotHistogram(allHistogramsCell, ...
              allTimeStampCell, ...
              histBinEdges, ...
              threshold, ...
              eventTimesStr, ...
              eventLabels);
