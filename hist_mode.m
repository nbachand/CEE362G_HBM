function [mode] = hist_mode(data,title_str)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figure()
h = histogram(data, 'Normalization', 'pdf');
title(title_str)

[~, max_id] = max(h.Values);
mode = mean(h.BinEdges(max_id:max_id+1));

end