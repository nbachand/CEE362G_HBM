%% Load data:
% xg : location along core
% y  : measured concentration

Cdata = dlmread( 'Cdata.dat' ) ;
xg = Cdata(:,1); y = Cdata(:,2); clear Cdata
%% Create H matrix
m = 19; n = 19;
H = zeros(n,n);
for i=1:19
    s = zeros(1,19);
    s(i) = 1;
    H(:,i) = Aquiclude_Concentration(s, 4.5, 0, 0);
end
%% Use MCMC to numerically compute the posterior density function
R0 = eye(19);
Q0 = exp(-(xg-xg').^2);
X = ones(m,1);

%Program the MCMC
PHI = H*X;
T = null(PHI')';
z = T*y;

C = @(th) T*(10^th(1)*H*Q0*H' + 10^th(2)*R0)*T';
L = @(th) -0.5*logdet(C(th)) -  0.5*z'*(C(th)\z);

N = 1000000; %length of the chain

%%
th_MC = zeros(N,2); %declaration/allocation with zeros
th_MC(1,:) = [5,2]; %initialization 
for k = 2:N
    %generate candidates
    th_MC(k,:) = th_MC(k-1,:) + [0.01*randn,0.01*randn];
    u = rand;
    if u>exp(L(th_MC(k,:))-L(th_MC(k-1,:)))
        th_MC(k,:) = th_MC(k-1,:);
    end
end

figure(1)
histogram(th_MC(:,1), 'Normalization', 'pdf')
title('\theta_1')
figure(2)
histogram(th_MC(:,2), 'Normalization', 'pdf')
title('\theta_2')

%%
% 
%  Numerically, the most likely theta1 is 5.265 and theta2 is 1.523 by
%  looking at the peaks, which is pretty much what we get from the previous
%  method! Below I plot the results and comparing with the first method, we
%  get basically the same results.
% 

%% Compare results from the two methods of inferring the hyperparameters
Q = 10^5.265* Q0;
R = 10^1.523* R0;
[s_est, SIG, LAMBDA, MU] = GenLinInv(y, H, R, X, Q);
stderr = sqrt(max(diag(SIG),0));
figure, plot(xg, s_est, xg, s_est+2*stderr, ':', xg, s_est-2*stderr, ':');
legend('Estimate', 'Upper bound', 'Lower bound');
title('Estimate with the credible intervals')

yrep = Aquiclude_Concentration(s_est, 4.5 , 0, 0)'
figure, plot(xg, y,'d', xg, yrep, 'b')
legend('Data', 'Reproduced')
title('Data Comparison')
