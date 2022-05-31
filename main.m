%% Assignment 2
% Nicholas Bachand

clear
clf
close all
%% Problem 1

data = importdata("data_snod_2021.txt");
% figure(); plot(data(:,1), data(:,2))

x = 1000;
u = 1;
D = 4E-2;
f = @(T) x./(2.*sqrt(pi.*D.*T.^3)).*exp(-((x-u.*T).^2)/(4.*D.*T));

x_test = 1:1:2000;
%figure(); plot(x_test, arrayfun(f, x_test))

%% 1)
meas_t = 1026:1:1249;
del_tau = 1;
signal_width = 300;
y = data(meas_t-data(1,1),2);

%figure(); plot(meas_t,y)

source_tau = 0:del_tau:meas_t(end)-x/u+signal_width/2;

fprintf('\n time window spans days %i to %i \n', source_tau(1), source_tau(end));
fprintf('\n delta tau is %i day \n', del_tau);

num_data = length(meas_t);
num_time_steps = length(source_tau);

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

H(1:10,1:10);

%% 2)
%% Q and R matrices
% Introduce a covariance matrix
theta1 = 1; theta2 = -8;
xg = source_tau;

q =@(x) -abs(x) ;
Q0 = toeplitz(q(xg));
R0 = eye(num_data);

Q0(1:10,1:10)
R0(1:10)

%% Metropolis-Hastings
%Program the 

X_mu = ones(num_time_steps,1);
th = [theta1, theta2];

theta = th_MC(end,:);

sprintf('Theta1: %g., Theta2: %g', round(theta(1),2), round(theta(2),2))

%% Updated Q and R

Q = (10^theta(1))*Q0;
R = (10^theta(2))*R0;

%% Inversion
[s,Sigma,LAMBDA,MU] = GenLinInv(y,H,R,X_mu,Q);
% figure()
% plot(xg,s)
%% Mean and Confidence Intervals
s_upp = s + 2*sqrt(diag(Sigma));
s_lower = s - 2*sqrt(diag(Sigma));

figure()
plot(xg,s,xg,s_lower,xg,s_upp)
title('mean and lower/upper confidence intervals')
ylabel("s")
xlabel("days")

%% 3)
ds = diff(s);
figure(); plot(xg(1:end-1), ds)

title("derivative of estimated s")

%% Hypothetical
hyp_true = zeros(num_time_steps,1);
hyp_true(50:70) = 1;

y_hyp = H*hyp_true;

[s_hyp,Sigma,LAMBDA,MU] = GenLinInv(y_hyp,H,R,X_mu,Q);

figure(); plot(xg, hyp_true, xg, s_hyp)

legend("true", "estimate")
title("hypothetical")
ylabel("s")
xlabel("days")




%% 4)
%% Hierarchical Log-Linear with Prior Information
X_mu = zeros(num_time_steps,3);
X_mu(:,1) = 1;
X_mu(62:82,2) = 1;
X_mu(112:122,3) = 1;

crit = 0.001;
maxstep = .1;
th_init = [theta1, theta2];
[theta,V,L] = HLLwPI(Q0,R0,X_mu,H,y,th_init,crit,maxstep);
sprintf('Theta1: %g., Theta2: %g', round(theta(1),2), round(theta(2),2))

%% Updated Q and R

Q = (10^theta(1))*Q0;
R = (10^theta(2))*R0;

%% Inversion
[s,Sigma,LAMBDA,MU] = GenLinInv(y,H,R,X_mu,Q);
% figure()
% plot(xg,s)
%% Mean and Confidence Intervals
s_upp = s + 2*sqrt(diag(Sigma));
s_lower = s - 2*sqrt(diag(Sigma));

figure()
plot(xg,s,xg,s_lower,xg,s_upp)
title('mean and lower/upper confidence intervals')
ylabel("s")
xlabel("days")

