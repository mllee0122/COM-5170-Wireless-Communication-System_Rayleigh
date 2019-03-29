%Wireless Communication Network
%Homework3_20181126
%107064522

clc;
clear;

%% 1. 
%  Rayleigh fading channel simulator based on the Filtered Gaussian Noise method
fmT = [0.01 0.1 0.5];
fm = fmT; % Assume simulation step size T = 1
omega_p = 1;    % Assume average power omega_p = 1
zeta = 2-cos(pi*fmT/2)-sqrt((2-cos(pi*fmT/2)).^2-1);
sigma = sqrt(((1+zeta)*omega_p)./((1-zeta)*2));
fm_tau = repmat((0:0.01:10), length(fmT), 1);
tau = fm_tau./(fm)';
t = linspace(0, 300, length(fm_tau));


w_1 = zeros(length(sigma), length(fm_tau)-1);
w_2 = zeros(length(sigma), length(fm_tau)-1);
for i = 1:length(sigma)
    w_1(i, :) = normrnd(0, sigma(i), 1, length(fm_tau)-1);  %  zero mean for Rayleigh fading
    w_2(i, :) = normrnd(0, sigma(i), 1, length(fm_tau)-1);  %  zero mean for Rayleigh fading
end

g_I = zeros(length(fm), length(fm_tau));
g_Q = zeros(length(fm), length(fm_tau));
for i = 1:length(fm_tau)-1
    g_I(:, i+1) = (zeta'.*g_I(:, i)) + ((1-zeta)'.*w_1(:, i));
    g_Q(:, i+1) = (zeta'.*g_Q(:, i)) + ((1-zeta)'.*w_2(:, i));
end

j = sqrt(-1);
g = g_I + j*g_Q;
r = abs(g); % Amplitude

% Autocorrelation (Calculation)
phi1 = zeros(length(fm), length(tau));
for i = 1:length(tau)
    phi1(:, i) = (((1-zeta)./(1+zeta)).*(sigma.^2))'.*(zeta'.^(tau(:, i)));
end

% Autocorrelation (Call Function)
phi2 = [autocorr(g(1,:), length(fm_tau)-1);
    autocorr(g(2,:), length(fm_tau)-1);
    autocorr(g(3,:), length(fm_tau)-1)];

%% Plot
figure(1)
plot(t, 10*log10(r))	% dB scale
title('Gaussian Noise Method')
xlabel('Time, t/T')
ylabel('Envelope Level (dB)')
legend({'fmT = 0.01', 'fmT = 0.1', 'fmT = 0.5'}, 'Location', 'southeast')

figure(2)
plot(t, phi1)
title('Gaussian Noise Method (Calculation)')
xlabel('Time Delay, fm\tau')
ylabel('Autocorrelation')
legend({'fmT = 0.01', 'fmT = 0.1', 'fmT = 0.5'}, 'Location', 'northwest')

figure(3)
plot(t, phi2)
title('Gaussian Noise Method (Call Function)')
xlabel('Time Delay, fm\tau')
ylabel('Autocorrelation')
legend({'fmT = 0.01', 'fmT = 0.1', 'fmT = 0.5'}, 'Location', 'northwest')