clear all
close all
clc

fontsize = 12.5;

Lag = [1.1 0.8 -0.5 -0.2 -0.8 -3.6 -3.9 -4.1 -4.5 -6.8 -11.1];
Lagstd = [0.7 0.8 0.8 0.7 0.8 0.9 0.8 1 1.1 1.2 1.4];
ATLs = [3 6 10 15 21 28 36 45 55 66 78];

LeLa_fit = polyfit(ATLs,Lag,1);
LeLa_fit = polyval(LeLa_fit,ATLs);

subplot(1,3,1)
errorbar(ATLs,Lag,Lagstd,'o-k','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'LineStyle','none')
hold on
plot(ATLs,LeLa_fit,'--','color',[0.5 0.5 0.5],'LineWidth',2)
grid on
ylim([-13 4])
xlim([-2 82])
grid on
set(gca,'FontSize',fontsize)
ylabel('Lag in Percentage of a 90 bmp beat')
xlabel('Transmission Latency (ms)')

% time paramters
fs = 10000;
T = 1/fs;
dur = 10;
t = 0:T:(dur-T);
ntime = length(t);

% metronome parameters
IOI = 0.667;
ft = 1/IOI;
bpm = ft*60;
x = exp(1i*2*pi*t*ft);

c1cycs = [0:2:99];
c1cycs = c1cycs(:);
c2cycs = [2:2:99];
c2cycs = c2cycs(:);
ATLs = [0 3 6 10 15 21 28 36 45 55 66 78]*fs/1000;
duet_start_n = 0;
idem = 1;
lets = {'B','C','D'};
a = 1;
b = -1;
A = 0.5;
F = 1;
D = 0.05;
LeLas = [];
f_delt = 0.02;

for ATL = ATLs; % the different ATLs to be tested
    
    zdel = round(0.222*fs); % delay as a fraction of the sampling rate
    
    % oscillator initialization
    zm1 = exp(1i*2*pi)*ones(size(t));
    f1 = ft;
    zm2 = exp(1i*2*pi)*ones(size(t));
    f2 = ft;
    
    c1on = zeros(size(t));
    c2on = zeros(size(t));
    
    c1 = 0;
    c2 = 0; % compensating for the initial response
    n = 1;
    n_switch = [1];
    c1_on = 1;
    
    while n<=zdel
        if n>ATL + 1
            zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - A*zm1(n));
            zm2(n+1) = zm2(n) + (f2 + f_delt)*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + F*zm1(n-ATL));
            n = n + 1;
        else
            zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - zm1(n)/abs(zm1(n)));
            zm2(n+1) = zm2(n) + (f2 + f_delt)*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2));
            n = n + 1;
        end
    end
    
    while n < ntime
        
        if (real(zm2(n) - zm2(n-1))) < 0 && (real(zm2(n-1) - zm2(n-2))) > 0 && n>(n_switch(end)+1000) && c1_on == 1 % true everytime a cycle ends
            c1_on = 0;
            n_switch = [n_switch n];
        end
        if (real(zm1(n) - zm1(n-1))) < 0 && (real(zm1(n-1) - zm1(n-2))) > 0 && n>(n_switch(end)+1000) && c1_on == 0 % true everytime a cycle ends
            c1_on = 1;
            n_switch = [n_switch n];
        end
        if c1_on
            c1on(n) = 1;
        end
        if ~c1_on
            c2on(n) = 1;
        end
        
        zm1(n+1) = zm1(n) + (f1 + c2on(n-ATL)*f_delt)*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - A*c1on(n)*zm1(n) + F*c2on(n-ATL)*zm2(n-ATL) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + (f2 + c1on(n-ATL)*f_delt)*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) - A*c2on(n)*zm2(n) + F*c1on(n-ATL)*zm1(n-ATL) - D*(1/f2).*(zm2(n - zdel)));
        
        n = n + 1;
    end
    
    [pks, rawbeats1] = findpeaks(real(zm1));
%         plot(real(zm1))
%         hold on
%         plot(rawbeats1,pks,'*')
    [pks, rawbeats2] = findpeaks(real(zm2));        
%         plot(real(zm2))
%         plot(rawbeats2,pks,'*')
%         plot(c1on)
%         length(rawbeats1)
%         length(rawbeats2)
%         pause
%         clf
    asynchronies = [];
    
    for iswitch=2:length(n_switch)
        if mod(iswitch,2) == 0
            all_asynchs = n_switch(iswitch) - rawbeats1;
            clipped_asynchs = min(abs(all_asynchs),fs/10);% clip asynchronies with magnitude greater than 0.1sec            
            ROI_asynchs = rawbeats1(find(clipped_asynchs<(fs/10)));
            asynchronies = [asynchronies n_switch(iswitch) - ROI_asynchs(1)];
        else
            all_asynchs = n_switch(iswitch) - rawbeats2;
            clipped_asynchs = min(abs(all_asynchs),fs/10);% clip asynchronies with magnitude greater than 0.1sec            
            ROI_asynchs = rawbeats2(find(clipped_asynchs<(fs/10)));
            asynchronies = [asynchronies n_switch(iswitch) - ROI_asynchs(1)];
        end
    end
    
    asynchs = asynchronies(1:2:end-1) + asynchronies(2:2:end);
    
    LeLas = [LeLas mean(asynchs)/10]
    
end
LeLa_all  = LeLas;

subplot(1,3,2)
plot(ATLs/10,-100*LeLa_all/666,'ok','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'LineStyle','none')
hold on
plot([3 6 10 15 21 28 36 45 55 66 78],LeLa_fit,'--','color',[0.5 0.5 0.5],'LineWidth',2)
ylim([-13 4])
xlim([-2 82])
grid on
set(gca,'FontSize',fontsize)
xlabel('Transmission Latency (ms)')

D = 0.36;
LeLas = [];

for ATL = ATLs; % the different ATLs to be tested
    
    zdel = round(0.222*fs); % delay as a fraction of the sampling rate
    
    % oscillator initialization
    zm1 = exp(1i*2*pi)*ones(size(t));
    f1 = ft;
    zm2 = exp(1i*2*pi)*ones(size(t));
    f2 = ft;
    
    c1on = zeros(size(t));
    c2on = zeros(size(t));
    
    c1 = 0;
    c2 = 0; % compensating for the initial response
    n = 1;
    n_switch = [1];
    c1_on = 1;
    
    while n<=zdel
        if n>ATL + 1
            zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - A*zm1(n));
            zm2(n+1) = zm2(n) + (f2 + f_delt)*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + F*zm1(n-ATL));
            n = n + 1;
        else
            zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - zm1(n)/abs(zm1(n)));
            zm2(n+1) = zm2(n) + (f2 + f_delt)*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2));
            n = n + 1;
        end
    end
    
    while n < ntime
        
        if (real(zm2(n) - zm2(n-1))) < 0 && (real(zm2(n-1) - zm2(n-2))) > 0 && n>(n_switch(end)+1000) && c1_on == 1 % true everytime a cycle ends
            c1_on = 0;
            n_switch = [n_switch n];
        end
        if (real(zm1(n) - zm1(n-1))) < 0 && (real(zm1(n-1) - zm1(n-2))) > 0 && n>(n_switch(end)+1000) && c1_on == 0 % true everytime a cycle ends
            c1_on = 1;
            n_switch = [n_switch n];
        end
        if c1_on
            c1on(n) = 1;
        end
        if ~c1_on
            c2on(n) = 1;
        end
        
        zm1(n+1) = zm1(n) + (f1 + c2on(n-ATL)*f_delt)*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - A*c1on(n)*zm1(n) + F*c2on(n-ATL)*zm2(n-ATL) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + (f2 + c1on(n-ATL)*f_delt)*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) - A*c2on(n)*zm2(n) + F*c1on(n-ATL)*zm1(n-ATL) - D*(1/f2).*(zm2(n - zdel)));
        
        n = n + 1;
    end
    
    [pks, rawbeats1] = findpeaks(real(zm1));
%         plot(real(zm1))
%         hold on
%         plot(rawbeats1,pks,'*')
    [pks, rawbeats2] = findpeaks(real(zm2));        
%         plot(real(zm2))
%         plot(rawbeats2,pks,'*')
%         plot(c1on)
%         length(rawbeats1)
%         length(rawbeats2)
%         pause
%         clf
    asynchronies = [];
    
    for iswitch=2:length(n_switch)
        if mod(iswitch,2) == 0
            all_asynchs = n_switch(iswitch) - rawbeats1;
            clipped_asynchs = min(abs(all_asynchs),fs/10);% clip asynchronies with magnitude greater than 0.1sec            
            ROI_asynchs = rawbeats1(find(clipped_asynchs<(fs/10)));
            asynchronies = [asynchronies n_switch(iswitch) - ROI_asynchs(1)];
        else
            all_asynchs = n_switch(iswitch) - rawbeats2;
            clipped_asynchs = min(abs(all_asynchs),fs/10);% clip asynchronies with magnitude greater than 0.1sec            
            ROI_asynchs = rawbeats2(find(clipped_asynchs<(fs/10)));
            asynchronies = [asynchronies n_switch(iswitch) - ROI_asynchs(1)];
        end
    end
    
    asynchs = asynchronies(1:2:end-1) + asynchronies(2:2:end);
    
    LeLas = [LeLas mean(asynchs)/10]
    
end
LeLa_all  = LeLas;

subplot(1,3,3)
plot(ATLs/10,-100*LeLa_all/666,'ok','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'LineStyle','none')
hold on
plot([3 6 10 15 21 28 36 45 55 66 78],LeLa_fit,'--','color',[0.5 0.5 0.5],'LineWidth',2)
ylim([-13 4])
xlim([-2 82])
grid on
set(gca,'FontSize',fontsize)
xlabel('Transmission Latency (ms)')
