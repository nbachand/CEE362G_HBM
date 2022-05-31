function [th_M] = untitled(H, X, Q0, R0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
PHI = H*X;
T = null(PHI')';
z = T*y;

C = @(th) T*(10^th(1)*H*Q0*H' + ...
    10^th(2)*R0)*T';
%Logarithm of the posterior PDF
L = @(th) -0.5*logdet(C(th)) - ...
    0.5*z'*(C(th)\z);

N = 100; %length of the chain
th_MC = zeros(N,2); %declaration/allocation with zeros
th_MC(1,:) = [1.1,-4.2]; %initialization 
for k = 2:N
    %generate candidates
    th_MC(k,:) = th_MC(k-1,:) + [0.1*randn,0.01*randn];
    u = rand;
    if u>exp(L(th_MC(k,:))-L(th_MC(k-1,:)))
        th_MC(k,:) = th_MC(k-1,:);
    end
end
end