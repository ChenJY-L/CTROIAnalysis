function [rows, cols] = detectEpidermis(currentFrame, edgeRoiMask, thresh)

edges = edge(currentFrame, 'canny', thresh);
edges = edges .* edgeRoiMask;
[edgeRowIdx, edgeColIdx] = find(edges == 1);
[~, ia] = unique(edgeColIdx, 'first');  % 每列第一个1的位置
rows = edgeRowIdx(ia);
cols = edgeColIdx(ia);
end