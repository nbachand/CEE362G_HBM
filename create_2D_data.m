clear, close all

% 3 sources at 3 locations s1, s2, s3. They each have mag/time release
% vector rel1, rel2, rel3 that tells us the ground truth. 
% we can calc concentration map for x,y grid given s_i and rel_i
% The concentration at each detector point is the sum of the 3
% concentration maps. -- ASSUMPTION. Further work would involve interplay
% of multiple sources.
%% Data
u = [0,0]; %Vel [m/day]
D = 1.9094; %Diffusion [m^2/day]
%impulse response function

f = @(s_l, d_l, T) diffusion_f(s_l, d_l, T, D, u);
% f = @(s_l, d_l, T) (norm(s_l-d_l))./sqrt(4*pi*D*T.^3).*exp(-(norm(d_l-s_l).^2./(4*D*T)));
% f = @(s_l, d_l, T) (norm(s_l-d_l))./sqrt(4*pi*D*T.^3).*exp(-(norm(d_l-s_l-u*T).^2./(4*D*T)));
% f = @(s_l, d_l, T) (norm(s_l-d_l))./sqrt(4.*pi.*D.*T.^3).*exp(-(vecnorm((d_l'-s_l').*ones(size(u'*T))-u'*T).^2./(4*D*T)));
%% Simulate data
D_loc = [5,8; 8,12; 12,6; 16,3; 14,14;10,2;20,1;3,5;1,15;18,10]; % detector locations
D_loc = [5,8; 8,12; 12,6];
S_loc = [4,10; 6,3; 10,12; 16,7;]; % source locations

n_sources = 4;
n_detectors = length(D_loc);

figure(1); hold on
indices = [5,8,9];
% scatter(D_loc(indices,1),D_loc(indices,2) );
scatter(D_loc(:,1),D_loc(:,2) );
scatter(S_loc(:,1),S_loc(:,2) );
xlim([0 20])
ylim([0 15])
xlabel('(m)'), ylabel('(m)');
legend('Detectors', 'Sources');
title('Detector and Source Locations');
%% times and values of actual release % times in days
tend = 200;
taus = (1:1:tend)'; [m,~] = size(taus);
background = 5;
s_gt = background*ones(m, n_sources); % 1 column per source

% True Release
for i=1:n_sources
    if i==1
        s_gt(61:80,i) = background+47*ones(20,1);
        s_gt(111:120,i) = background+19*ones(10,1);
    elseif i==2
        small = [20*ones(1,20), 5*ones(1,20)];
        rep = repmat(small,1,round(tend/length(small)));
        s_gt(:,i) = rep(1:tend);
    elseif i==3
%         s_gt(41:80,i) = background+35*ones(40,1);
    elseif i==4
        s_gt(:,i) = s_gt(:,i)*2;
    end
end

%the release analytically
san1 = @(t) background + 47*(t>60).*(t<=80) + 19*(t>110).*(t<=120);
san2 = @(t) 20*(t>0).*(t<=20)+5*(t>20).*(t<=40) + 20*(t>40).*(t<=60)+5*(t>60).*(t<=80) + 20*(t>80).*(t<=100) +5*(t>100).*(t<=120)+ 20*(t>120).*(t<=140) +5*(t>140).*(t<=160)+ 20*(t>160).*(t<=180) +5*(t>180).*(t<=200);
% san3 = @(t) background+ 35*(t>40).*(t<=80);
% san2 = @(t) background*ones(size(t));
san3 = @(t) background*ones(size(t));
san4 = @(t) background*2*ones(size(t));

figure(2)

plot(taus,s_gt(:,1),taus,s_gt(:,2),taus,s_gt(:,3), 'm', taus,s_gt(:,4), 'c', taus,san1(taus),'b:', taus,san2(taus), 'r:', taus,san3(taus), 'm:', taus,san3(taus), 'c:')
title('"True" release')
xlabel('time (days)')
ylabel('release')
legend('source 1', 'source 2', 'source 3','source 4','source 1 analytical', 'source 2 analytical', 'source 3 analytical', 'source 4 analytical');
% legend('source 1','source 1 analytical');

%% concentration observations
n = 221;    %number of observations
ts = 210-200;  %start time
te = 430-200;  %end time
t = (ts:1:te)'; %times of sampling
s1 = S_loc(1,:); s2 = S_loc(2,:); s3 = S_loc(3,:); s4 = S_loc(4,:); 
conc = zeros(n, n_detectors); %per detector, we get samples at different times

% tau is the source times
% tau = 10:.01:400;
for niter = 1:n % for each time sample
    for i=1:n_detectors % for each of the n_detectors, find the concentrations
        F1 = @(tau) f(s1, D_loc(i,:),t(niter)-tau).*san1(tau);
        F2 = @(tau) f(s2, D_loc(i,:),t(niter)-tau).*san2(tau);
        F3 = @(tau) f(s3, D_loc(i,:),t(niter)-tau).*san3(tau);
        F4 = @(tau) f(s4, D_loc(i,:),t(niter)-tau).*san4(tau);
        conc(niter,i) = quad(F1,0.5,200.5) + quad(F2,0.5,200.5)+ quad(F3,0.5,200.5)+quad(F4,0.5,200.5);
    end
end

rng(1);
y = conc+0.2*randn(n,n_detectors);
%%
figure(3) 
hold on

for i=1:n_detectors
    plot(t, conc(:,i), t, y(:,i) ,'.');
end

title('Detector Measurements');
legend(repmat([{'true signal'},{'sampled'}], 1, n_detectors));
xlabel('time (days)');
hold off;

%% SAVE
methane_data = [t,y];
save methane_data.txt -ascii methane_data
% end simulate data 

%% MISC
%% -------------------- create a video -------------------------------
% Concentration map from 1 s_i
% x = (0:0.01:20);
% y = (0:0.01:15);
% tsource = (1:5:100);
% z1 = zeros(length(x), length(y), length(tsource));
% s1 = [5,5];
% for i=1:length(x)
%     for j=1:length(y)
%         detector_loc = [x(i), y(j)];
%         for k=1:length(tsource)
%             z1(i,j,k) = f(s1, detector_loc, tsource(k));
%         end
%     end
% end
% z = z1;
% v = VideoWriter('test.avi');
% v.FrameRate = 5; 
% open(v);
%  
% set(gca,'nextplot','replacechildren', 'YDir','normal', 'visible', 'off');
% colorbar;
% 
% for k = 1:length(tsource) 
%    imagesc(z(:,:,k));  
%    caxis([0, 0.01]);
%    frame = getframe;
%    writeVideo(v,frame);
% end
% close(v);
% 

%-------------------------------------------------------