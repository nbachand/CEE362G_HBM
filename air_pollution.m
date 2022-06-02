%% Get estimates from data
clear, close all

u = [0,0]; %Vel [m/day] (of wind or something. We will try both 0 velocity so it's just diffusion and we will try constant velocity wind)
D = 1.9094;  %Diffusion[m^2/day] | % methane is  D = 0.221 cm^2/s or 1.9094 m^2/day
f = @(s_l, d_l, T) diffusion_f(s_l, d_l, T, D, u);

%% Load data
load methane_data.txt
ty = methane_data(:,1);
y = methane_data(:,2:end); 
[n_measurements, n_detectors] = size(y); 
n_sources = 4;
y = y(:); % stack columns of measurements on top of each other into single column

%% Select time period for source strength
tsource = 1:200; %We will try to find each source strength, s_i, at these times (in days)
m = length(tsource);

%% Select detector and source locations
D_loc = [5,8; 8,12; 12,6; 16,3; 14,14;10,2;20,1;3,5;1,15;18,10]; % detector locations
S_loc = [4,10; 6,3; 10,12; 16,7;];

figure(1); hold on
scatter(D_loc(:,1),D_loc(:,2) );
scatter(S_loc(:,1),S_loc(:,2) );
legend('Detectors', 'Sources');
title('Detector and Source Locations');
xlabel('(m)'), ylabel('(m)');
xlim([0 20])
ylim([0 15])
hold off
figure(1); hold on
scatter(S_loc(:,1),S_loc(:,2), 'r' );
title('Source Locations');
xlabel('(m)'), ylabel('(m)');
xlim([0 20])
ylim([0 15])
hold off


%% Construct H matrix for single detector
H = zeros(n_detectors*length(ty), n_sources*length(tsource));
for di=1:n_detectors % for each detector, create sub H matrix
    Hi = zeros(n_measurements, n_sources*length(tsource)); 
    for si=1:n_sources
        source_loc = S_loc(si,:);
        det_loc = D_loc(di,:);
        Hi(:,(si-1)*m+1:(si)*m) = detector_H(ty, tsource, source_loc, det_loc, u, D);
    end
    H((di-1)*n_measurements+1:di*n_measurements,:) = Hi;  
end
hal_H = H;

[U,S,V] = svd(H,'econ');  
sigma = diag(S); nrH = sum(sigma/sigma(1)>eps); %numerical rank
 
%% Start inversion
[n,m] = size(H);
Q0 = toeplitz(-(0:1:m-1));
X = [ones(m,1)]; 
R0 = eye(n);

[thetaLL,VLL,L] = HLLwPI(Q0,R0,X,H,y,[0;-1],0.001,0.5);
th = 10.^thetaLL;
Q = th(1)*Q0;
R = th(2)*R0;
[s_hat,V,LAMBDA,MU] = GenLinInv(y,H,R,X,Q);
%% Plot reconstructed signals
figure(4)
sig = sqrt(diag(V));

t_n  = length(tsource);
hold on;
subplot(2,2,1);
plot(tsource,s_hat(1:t_n,1), tsource,s_hat(1:t_n,1)+2*sig(1:t_n,1),'--',tsource,s_hat(1:t_n,1)-2*sig(1:t_n,1),'--' );
legend('Reconstructed','Lower Bound','Upper Bound');
title('Source 1 Reconstruction');
ylim([-1 100]);

subplot(2,2,2);
plot(tsource,s_hat(t_n+1:2*t_n,1), tsource,s_hat(t_n+1:2*t_n,1)+2*sig(t_n+1:2*t_n,1),'--',tsource,s_hat(t_n+1:2*t_n,1)-2*sig(t_n+1:2*t_n,1),'--' );
legend('Reconstructed','Lower Bound','Upper Bound');
title('Source 2 Reconstruction');
ylim([-1 100]);

subplot(2,2,3);
plot(tsource,s_hat(2*t_n+1:3*t_n,1), tsource,s_hat(2*t_n+1:3*t_n,1)+2*sig(2*t_n+1:3*t_n,1),'--',tsource,s_hat(2*t_n+1:3*t_n,1)-2*sig(2*t_n+1:3*t_n,1),'--' );
legend('Reconstructed','Lower Bound','Upper Bound');
title('Source 3 Reconstruction');
ylim([-1 100]);
hold off;

subplot(2,2,4);
plot(tsource,s_hat(3*t_n+1:4*t_n,1), tsource,s_hat(3*t_n+1:4*t_n,1)+2*sig(3*t_n+1:4*t_n,1),'--',tsource,s_hat(3*t_n+1:4*t_n,1)-2*sig(3*t_n+1:4*t_n,1),'--' );
legend('Reconstructed','Lower Bound','Upper Bound');
title('Source 4 Reconstruction');
ylim([-1 100]);
hold off;
