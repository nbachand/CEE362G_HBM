function [H] = satellite_H(true_s, kernel_size)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n_pixels = numel(true_s);
n_data = n_pixels/kernel_size;

% initiate H
H = zeros(n_data,n_pixels);

% fill in h to act as a 1D average pooling layer
for i = 1:kernel_size
    H(:,i:i+n_data) = H(:,i:i+n_data)+eye(n_data)./kernel_size;
end