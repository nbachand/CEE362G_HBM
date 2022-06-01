% Air Pollution
clear, close all

sanity_check = true;
% 2D
u = [0,0]; %Vel [m/day] (of wind or something. We will try both 0 velocity so it's just diffusion and we will try constant velocity wind)
D = 1.9094;  %Diffusion[m^2/day] | % methane is  D = 0.221 cm^2/s or 1.9094 m^2/day

if sanity_check
    % 1D -- sanity check
    u = [1,0]; %Vel [m/day]
    D = 4.E-2; %Diffusion [m^2/day]    
end


% impulse response function
% s_l is a source location (x,y) on grid 
% d_l is a detector location (x,y) on grid
f = @(s_l, d_l, T) (norm(s_l-d_l))./sqrt(4*pi*D*T.^3).*exp(-(norm(d_l-s_l-u*T).^2./(4*D*T)));

tau = (1:1:500)';
out_tau = zeros(size(tau));
for i=1:length(tau)
    out_tau(i) = f([0,0],[40,0],tau(i));
end
% plot(tau, out_tau)

%% Load data
load data_snod_2021.txt
ty = data_snod_2021(:,1);
% measurements at different detectors at different times
y = data_snod_2021(:,2:end); 
[n_measurements, n_detectors] = size(y); % should be 10 detectors and 3 sources
n_sources = 1;

y = y(:); % stack columns of measurements on top of each other into single column

%% Select time period for source strength
tsource = 1:270; %We will try to find each source strength, s_i, at these times (in days)
m = length(tsource);

%% Select detector and source locations
D_loc = [10,2; 3,5; 20,1; 1,15; 5,8; 14,14; 16,3; 8,12; 12,6; 18,10];
S_loc = [6,3; 10,12; 16,7];

if sanity_check
    D_loc = [1000,0];
    S_loc = [0,0];
end


figure(1); hold on
scatter(D_loc(:,1),D_loc(:,2) );
scatter(S_loc(:,1),S_loc(:,2) );
legend('Detectors', 'Sources');
title('Detector and Source Locations');
xlabel('(m)'), ylabel('(m)');
hold off

%% Construct H matrix for single detector
% s_est = zeros(n_sources*length(tsource),1); % stack S1 (t0 - tn) on S2 (t0 -tn) on S3 (t0-tn)

H = zeros(n_detectors*length(ty), n_sources*length(tsource));
for di=1:n_detectors % for each detector, create sub H matrix
    Hi = zeros(n_measurements, n_sources*length(tsource)); 
    for si=1:n_sources
        source_loc = S_loc(si,:)
        det_loc = D_loc(di,:)
        u, D
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

[thetaLL,VLL,L] = HLLwPI(Q0,R0,X,H,y,[0;-1],0.001,0.5)
th = 10.^thetaLL
Q = th(1)*Q0;
R = th(2)*R0;

[s_hat,V,LAMBDA,MU] = GenLinInv(y,H,R,X,Q);

figure(4)
sig = sqrt(diag(V));
plot(tsource,s_hat,tsource,s_hat+2*sig,'--',tsource,s_hat-2*sig,'--')

%% Analysis of results

figure(5)
ds_hat = diff(s_hat);
plot(1:m-1,ds_hat)
title('Slope of estimate')

%visually identify large-value pairs 60 to 80, 110 to 122
a = zeros(m,1);
a(60:80)= 1; 
b = zeros(m,1);
b(110:122) = 1;
% Explore whether this a and b could be resolved with previous model
A = [ones(m,1),Q*H'];
W = orth(A);
norm(W'*a)^2/norm(a)^2
norm(W'*b)^2/norm(b)^2
[a_hat,~,~,~] = GenLinInv(H*a,H,R,X,Q);
norm(a_hat)/norm(a)
figure(6)
plot(1:m,a,1:m,a_hat)
legend('true','estimate')

[b_hat,~,~,~] = GenLinInv(H*b,H,R,X,Q);
norm(a_hat)/norm(a)
figure(7)
plot(1:m,b,1:m,b_hat)
legend('true','estimate')

%% A model with "zones"
X = [ones(m,1), a, b]; 
X = orth(X);

[thetaLL,VLL,L] = HLLwPI(Q0,R0,X,H,y,[0;-2],0.001,0.5)
th = 10.^thetaLL

Q = th(1)*Q0;
R = th(2)*R0;

[s_hat,V,LAMBDA,MU] = GenLinInv(y,H,R,X,Q);
figure(10)
sig = sqrt(diag(V));
plot(tsource,s_hat,tsource,s_hat+2*sig,'--',tsource,s_hat-2*sig,'--')

%end