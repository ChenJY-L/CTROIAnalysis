function [dicomFiles, filePaths] = scanDicomFiles(dicomDir)

fileList = dir(dicomDir);
dicomFiles = {};    % Store valid DICOM file names
filePaths = {};     % Store valid DICOM file paths

fprintf('正在扫描文件夹: %s\n', dicomDir);
for i = 1:length(fileList)
    fileName = fileList(i).name;
    filePath = fullfile(dicomDir, fileName);
    if fileList(i).isdir || startsWith(fileName, '.') || endsWith(fileName, [".mat", ".png", ".csv", ".avi"])
        continue;
    end
    try
        % Use dicominfo to quickly check if it's likely a DICOM
        % dicomInfo = dicominfo(filePath);
        dicomFiles{end+1} = fileName;
        filePaths{end+1} = filePath;
    catch
        % Not a DICOM file or error reading info, skip
        continue;
    end
end
end