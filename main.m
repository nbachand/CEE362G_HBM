%% Project
% Nicholas Bachand and Haley So

clear
clf
close all
%% Loading Data and Setting Parameters
data = importdata("data_snod_2021.txt");
%figure(); plot(data(:,1), data(:,2))

x = 1000;
u = 1;
D = 4E-2;
f = @(T) x./(2.*sqrt(pi.*D.*T.^3)).*exp(-((x-u.*T).^2)/(4.*D.*T));

x_test = 1:1:2000;
%figure(); plot(x_test, arrayfun(f, x_test))

%% Measurement Window and Frequency
meas_t = 1026:1:1249;
del_tau = 1;
signal_width = 300;
y = data(meas_t-data(1,1),2);

%figure(); plot(meas_t,y)

source_tau = 0:del_tau:meas_t(end)-x/u+signal_width/2;

fprintf('\n time window spans days %i to %i \n', source_tau(1), source_tau(end));
fprintf('\n delta tau is %i day \n', del_tau);

num_time_steps = length(source_tau);
num_data = length(meas_t);

%% Building advection-diffusion H

H = zeros(num_data, num_time_steps);

i = 1;
for t = meas_t
    j = 1;
    for tau = source_tau
        if t>tau
            H(i,j) = f(t-tau);
        end
        j = j+1;
    end
    i = i +1;
end

%% Adding on sattellite measurements (average pooling)

kernel_size = 4;
% added gaussian noise with stddev of 0.01. Output size is 56 to match Nick's H pooling matrix.
[H_pool, x_sampled, s_out] = average_H(meas_t, y, kernel_size, 56, 0.01);
% H_pool = satellite_H(y,4);
H = H_pool*H;
y = H_pool*y;


num_data = length(y);


%% Prior Assumptions (Q, R, X, theta
% Introduce a covariance matrix
theta1 = 1; theta2 = -1;
xg = source_tau;

q =@(x) -abs(x) ;
Q0 = toeplitz(q(xg));
R0 = eye(num_data);

Q0(1:10,1:10)
R0(1:10)

X_mu = zeros(num_time_steps,3);
X_mu(:,1) = 1;
X_mu(floor(62/kernel_size):ceil(82/kernel_size),2) = 1;
X_mu(floor(112/kernel_size):ceil(122/kernel_size),3) = 1;

th = [theta1, theta2];

%% MCMC Parameter Search Adv-Diff Only

N = 1000; %length of the chain
th_MC = MCMC(H, X_mu, Q0, R0, y, N, th);

%% Get Theta from Hist Mode
theta(1) = hist_mode(th_MC(:,1), '\theta_1');
theta(2) = hist_mode(th_MC(:,2), '\theta_1');
sprintf('Theta1: %g., Theta2: %g', round(theta(1),2), round(theta(2),2))


%% Updated Q and R
Q = (10^theta(1))*Q0;
R = (10^theta(2))*R0;

%% Inversion
[s,Sigma,LAMBDA,MU] = GenLinInv(y,H,R,X_mu,Q);

%% Mean and Confidence Intervals
s_upp = s + 2*sqrt(diag(Sigma));
s_lower = s - 2*sqrt(diag(Sigma));

figure()
plot(xg,s,xg,s_lower,xg,s_upp)
title('mean and lower/upper confidence intervals')
ylabel("s")
xlabel("days")

%% Plot ds
ds = diff(s);
figure(); plot(xg(1:end-1), ds)

title("derivative of estimated s")

