function roiMask = fixROI(oriROIMask, epidermisDiff)
[y, x] = find(oriROIMask);
y = y - epidermisDiff;
roiMask = false(size(oriROIMask));
idx = sub2ind(size(oriROIMask), y, x);
roiMask(idx) = true;

end