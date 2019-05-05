clear all
clc

figure('units','normalized','outerposition',[0 0 1 1])
fontsize = 15;

%%% Figure 1: Behavioral data from Repp & Doggett (2007) and %%%
%%% regression lines with 95% confidence interval.           %%%

IOIs = 1000:250:3500;

nmNMA = 1000*[-0.036 -0.048 -0.051 -0.066 -0.0775 -0.077 -0.103 -0.111 -0.1065 -0.142 -0.1475]; % From Repp 2007
nmNMAse = 1000*[0.008 0.012 0.013 0.015 0.016 0.015 0.013 0.021 0.022 0.025 0.030];
mNMA = 1000*[-0.0045 -0.008 -0.014 -0.014 -0.003 -0.007 0.000 -0.009 -0.006 -0.031 -0.031];
mNMAse = 1000*[0.0025 0.003 0.0055 0.008 0.008 0.007 0.011 0.0085 0.016 0.0125 0.016];

X = [ones(size(IOIs')) IOIs'];
y = mNMA';
[b,bint] = regress(y,X);
yhat = b(1)+b(2)*IOIs;

figure(1)
subplot(2,4,[1 2 5 6])
errorbar(IOIs,mNMA,mNMAse,'ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2,'DisplayName','Musicians')
hold on
errorbar(IOIs,nmNMA,nmNMAse,'ok','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',2,'DisplayName','Nonmusicians')
legend('show','Location','southwest');
plot(IOIs,yhat,'k-','linewidth',2);
ylim([-200 50])

X = [ones(size(IOIs')) IOIs'];
y = nmNMA';
[b,bint] = regress(y,X);
yhat = b(1)+b(2)*IOIs;

subplot(2,4,[1 2 5 6])
plot(IOIs,yhat,'k--','linewidth',1.5);
ylim([-200 50])
xlim([900 3600])
grid on
grid minor
xlabel('Stimulus Period (ms)','FontSize',fontsize)
ylabel('Asynchrony ± SE (ms)','FontSize',fontsize)
set(gca,'FontSize',fontsize)

% time paramters
fs = 10000;
T = 1/fs;
dur = 100; % seconds
t = 0:T:(dur+T);
ntime = length(t);

% stimulus periods to be tested
IOIs = 1.0:0.25:3.5;

% oscillator initialization and parameters
a = 1;
b = -1;
zdel = round(0.222*fs); % SMS delay as a fraction of the sampling rate
zm = 0.5*exp(1i*2*pi)*ones(size(t));
A = -0.5;

cols = {[0.0 0.6 0.2],[1.0 0.8 0.2]};
icolor = 1;

for D = [0.05 0.36];
    NMAs = [];
    NMAstds = [];
    for IOI = IOIs; % in seconds
        ft = 1/IOI;
        bps = ft;
        bpm = bps*60;
        x = exp(1i*2*pi*t*ft);
        F = 1;
        
        % initialize frequency vector
        f = ft;
        
        % main loop
        for n = 1:ntime-1
            
            if n > zdel
                
                zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)) - D*(1/(f)).*(zm(n - zdel)));
                
            else
                zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)));
            end
            
        end
        
        % finding the peaks and the NMA
        [pks, rawbeats1] = findpeaks(real(x));
        [pks, rawbeats2] = findpeaks(real(zm));
        rawbeats1 = [1 rawbeats1]; % adding the peak at n  = 1
        NMA = rawbeats2 - rawbeats1;
        NMAs = [NMAs mean(NMA(floor(end/2):end))];
        
    end
    subplot(2,4,[3 4 7 8])
    plot(1000*IOIs, 1000*NMAs/fs,'o','MarkerSize',10,'LineWidth',2,'DisplayName',sprintf('D=%.3f',D),'color',cols{icolor},'MarkerFaceColor',cols{icolor})        
    grid on
    grid minor
    hold on
    xlim([1000*IOIs(1) - 100 1000*IOIs(end) + 100])
    ylim([-200 50])
    xlabel('Stimulus Period (ms)','FontSize',fontsize)
    ylabel('Asynchrony (ms)','FontSize',fontsize)
    set(gca,'FontSize',fontsize)
    icolor = icolor + 1;
end

X = [ones(size(IOIs')) 1000.*IOIs'];
y = mNMA';
[b,bint] = regress(y,X);
yhatm = b(1)+b(2)*1000*IOIs;

X = [ones(size(IOIs')) 1000.*IOIs'];
y = nmNMA';
[b,bint] = regress(y,X);
yhatn = b(1)+b(2)*1000*IOIs;

plot(1000.*IOIs,yhatm,'k','linewidth',2,'DisplayName','Musicians');
plot(1000.*IOIs,yhatn,'k--','linewidth',2,'DisplayName','Nonmusicians');
legend('show','Location','southwest')
%%
% for i=2:3
%     subplot(3,4,i)
%     plot(IOIs,yhat,'k--','linewidth',2);
%     hold on
% end

