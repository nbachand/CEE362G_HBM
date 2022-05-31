function [H_sampled, x_sampled, s_out] = average_H(x, s, N, out_size, stddev)
% Computes H for averaging sample_size neighboring entries with or without
% noise and then downsamples to specified number out_size

% x: input locations
% s: values at input locations
% N: size of selection to average
% out_size: number of samples we want
% stddev: stddev of noise e.g. 0.01 (0 mean gaussian noise) 

n_in = numel(s);
H = zeros(n_in, n_in);

%---------- Create diagonal Matrix -------------
for i=1:N
    diagonal_matrix = diag(ones(1,n_in), i-1);
    H = H + diagonal_matrix(1:n_in, 1:n_in);
end
%    make sure each row's weights sum to 1
for row=1:n_in
    H(row,:) = H(row,:)/nnz(H(row,:));
end

%------------------Add noise--------------------
noise = stddev*randn(n_in);
H = H +noise;

%-----------------Downsample--------------------
to_sample = downsample((1:1:n_in), cast(n_in/out_size, 'uint8'))
x_sampled = x(to_sample);
H_sampled = H(to_sample,:);
x_sampled = x(1:out_size);
H_sampled = H_sampled(1:out_size,:);

s_out = H_sampled*s;