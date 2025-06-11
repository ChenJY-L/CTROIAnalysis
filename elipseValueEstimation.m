function [elipseCoeffs, elipseStat] = elipseValueEstimation(frame, roiMask, edges, unitY)
elipseValue = mean(frame .* uint8(roiMask), 1);
fitX = (find(elipseValue ~= 0) - mean(edges)).*unitY;
elipseValue = elipseValue(elipseValue ~= 0);
[elipseCoeffs, elipseStat] = robustfit(fitX, elipseValue);
end