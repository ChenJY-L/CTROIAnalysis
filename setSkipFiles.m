function filesToProcessIdx = setSkipFiles(numFiles)

promptSkip = {'请输入开头需要跳过的文件数量:', '请输入结尾需要跳过的文件数量:'};
dlgtitleSkip = '跳过文件';
definputSkip = {'0', '0'};
inputResponseSkip = inputdlg(promptSkip, dlgtitleSkip, [1 35], definputSkip);
if isempty(inputResponseSkip)
    disp('跳过文件数量输入取消。程序退出。');
    return;
end

numFilesToSkip(1) = str2double(inputResponseSkip{1});
numFilesToSkip(2) = str2double(inputResponseSkip{2});

fprintf('将跳过前 %d 个文件进行分析。\n', numFilesToSkip(1));
fprintf('将跳过后 %d 个文件进行分析。\n', numFilesToSkip(2));

% Determine the range of files to actually process
filesToProcessIdx = (numFilesToSkip(1) + 1) : (numFiles - numFilesToSkip(2));

end
