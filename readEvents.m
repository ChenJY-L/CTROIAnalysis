function [eventTimesStr, eventLabels] = readEvents(dicomDir)

filePath = fullfile(dicomDir, "events.csv");
data = readtable(filePath, 'Delimiter', ',');
eventTimesStr = cellstr(data.Time);
eventLabels = cellstr(data.Label);
end