%Wireless Communication Network
%Homework3_20181126
%107064522
 
clc;
clear;
 
%% 2. 
% Rayleigh fading channel simulator based on the Sum of Sinusoids method
fmT = [0.01 0.1 0.5];
fm = fmT; % Assume simulation step size T = 1
j = sqrt(-1);
 
alpha = 0;
g = zeros(2, 1001, length(fm));
for M = [8 16]
    N = 4*M+2;
    n = 1:M;
    beta = pi*n/M;
    fn = fm'*cos(2*pi*n/N);
    g_I = zeros(length(fm), length(n));
    g_Q = zeros(length(fm), length(n));
    gI = zeros(2, 1001, length(fm));
    gQ = zeros(2, 1001, length(fm));
    for t = 0:1000
        for i = 1:length(fm)
            g_I(i, :) = 2*cos(beta).*cos(2*pi*fn(i, :)*t);
            g_Q(i, :) = 2*sin(beta).*cos(2*pi*fn(i, :)*t);
        end
        if M == 8
            k = 1;
        elseif M == 16
            k = 2;
        end
        
    gI(k,t+1,:) = 2*sum(g_I,2)'+2^0.5*cos(alpha)*cos(2*pi*fm*t);
    gQ(k,t+1,:) = 2*sum(g_Q,2)'+2^0.5*sin(alpha)*cos(2*pi*fm*t);
    g(k,t+1,:) = sqrt(2)*cell2mat({complex(gI(k, t+1, :), gQ(k, t+1, :))});
    end
end
 
%% Plot
time = 0:300;
 
% Channel Output
for i = [1 2]
    figure(i)
    hold on
    for k = 1:length(fm)
        plot(time',10*log10(abs(g(i,time+1,k))))
    end
    
    if i == 1
        title('Sum of Sinusoids Method (M = 8)')
        xlabel('Time, t/T')
        ylabel('Envelope Level (dB)')
        legend({'fmT = 0.01', 'fmT = 0.1', 'fmT = 0.5'}, 'Location', 'southwest')
    else
        title('Sum of Sinusoids Method (M = 16)')
        xlabel('Time, t/T')
        ylabel('Envelope Level (dB)')
        legend({'fmT = 0.01', 'fmT = 0.1', 'fmT = 0.5'}, 'Location', 'northeast')
    end    
end
 
% Autocorrelation
for i = [1 2]
    figure(i+2)
    hold on
    for k = 1:length(fm)
        plot((0:fmT(k):10),autocorr(g(i,:,k),length((0:fmT(k):10))-1))
    end
    
    if i == 1
        title('Sum of Sinusoids Method (M = 8)')
        xlabel('Time Delay, fm\tau')
        ylabel('Autocorrelation')
        legend({'fmT = 0.01', 'fmT = 0.1', 'fmT = 0.5'}, 'Location', 'southwest')
    else
        title('Sum of Sinusoids Method (M = 16)')
        xlabel('Time Delay, fm\tau')
        ylabel('Autocorrelation')
        legend({'fmT = 0.01', 'fmT = 0.1', 'fmT = 0.5'}, 'Location', 'northeast')
    end    
end