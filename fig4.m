% Behavioral data from Chafe et al. 2010
% Lag as percentage of a 90 bpm beat

clear all
fontsize = 12.5;

% time paramters
fs = 1000;
T = 1/fs;
dur = 5.5;
t = 0:T:(dur-T);
ntime = length(t);

x = exp(1i*2*pi*t(1:2000));

z1 = ones(size(t));
for n=1:2000-1
    if n>200
        z1(n+1) = z1(n) + T*(z1(n)*(1 + 1i*2*pi - abs(z1(n))^2) + x(n) - z1(n-100));
    else
        z1(n+1) = z1(n) + T*(z1(n)*(1 + 1i*2*pi - abs(z1(n))^2) + x(n));
    end
end

z2 = ones(size(t));
for n=2000:3000-1
    
    z1(n+1) = z1(n) + T*(z1(n)*(1 + 1i*2*pi - abs(z1(n))^2));
    z2(n+1) = z2(n) + 0.9*T*(z2(n)*(1 + 1i*2*pi - abs(z2(n))^2));
    
end


for n=3000:4000-1
    
    z1(n+1) = z1(n) + 0.8*T*(z1(n)*(1 + 1i*2*pi - abs(z1(n))^2));
    z2(n+1) = z2(n) + 1.2*T*(z2(n)*(1 + 1i*2*pi - abs(z2(n))^2));
    
end

for n=4000:5000-1
    
    z1(n+1) = z1(n) + 1.2*T*(z1(n)*(1 + 1i*2*pi - abs(z1(n))^2));
    z2(n+1) = z2(n) + 0.8*T*(z2(n)*(1 + 1i*2*pi - abs(z2(n))^2));
    
end

for n=5000:5500-1
    
    z1(n+1) = z1(n) + 0.8*T*(z1(n)*(1 + 1i*2*pi - abs(z1(n))^2));
    z2(n+1) = z2(n) + 1.2*T*(z2(n)*(1 + 1i*2*pi - abs(z2(n))^2));
    
end

back = zeros(1,length(t),3);

back(:,2001:3000,3) = ones(1,1000,1);
back(:,3001:4000,1) = ones(1,1000,1);
back(:,4001:5025,3) = ones(1,1025,1);
back(:,5026:end,1) = ones(1,475,1);

[locs,beats1] = findpeaks(real(z1));
[locs,beats2] = findpeaks(real(z2(2000:end)));

time1 = linspace(2.5,6.5,length(z1));
time1 = time1(beats1);

time2 = linspace(3.95,6.5,length(z2(2000:end)));
time2 = time2(beats2);

subplot(2,2,1:2)
hold on
h = imagesc(back+0.85);
set(h, 'XData', [2.5, 6.5],'YData',[-2 2]);
plot(linspace(2.5,6.5,length(z1)),real(z1),'LineWidth',2,'color','b')
plot(linspace(3.95,6.5,length(z2(2000:end))),real(z2(2000:end)),'LineWidth',2,'color','r')
plot(linspace(2.5,3.95,length(x(1:2000))),real(x),'--','color',[0.5 0.5 0.5],'LineWidth',1)
grid on
tx = text(0.02,0.98,'B','Units', 'Normalized', 'VerticalAlignment', 'Top');
tx.FontSize = 20;
stem(time1,1.15*ones(size(time1)),'v','filled','LineStyle','none','MarkerSize',7,'color','b')
stem(time2,1.15*ones(size(time2)),'v','filled','LineStyle','none','MarkerSize',7,'color','r')
ylim([-1.2 1.4])
xlim([2.5 6.5])
set(gca,'FontSize',fontsize)
ylabel('Amplitude')
xlabel('time (s)')

Lag = [1.1 0.8 -0.5 -0.2 -0.8 -3.6 -3.9 -4.1 -4.5 -6.8 -11.1];
Lagstd = [0.7 0.8 0.8 0.7 0.8 0.9 0.8 1 1.1 1.2 1.4];
ATLs = [3 6 10 15 21 28 36 45 55 66 78];

LeLa_fit = polyfit(ATLs,Lag,1);
LeLa_fit = polyval(LeLa_fit,ATLs);

subplot(2,2,3)
errorbar(ATLs,Lag,Lagstd,'o-k','MarkerSize',8,'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1.25)
hold on
tx = text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top');
tx.FontSize = 20;
plot(ATLs,LeLa_fit,'color',[0.5 0.5 0.5],'LineWidth',2)
grid on
ylim([-13 4])
xlim([-2 82])
grid on
set(gca,'FontSize',fontsize)
ylabel('Lag in Percentage of a 90 bmp beat')
xlabel('Transmission Latency (ms)')
%%
%%% Lag as a function of ATL at the level of quarter notes

% time paramters
fs = 10000;
T = 1/fs;
dur = 5;
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
ATLs = [0 3 6 10 15 21 28 36 45 66 78]*fs/1000;
duet_start_n = 0;
idem = 1;
lets = {'B','C','D'};
a = 1;
b = -1;
A = 0.5;
F = 1;
D = 0.02;
LeLa_all = [];
for replicates = 1:10
    LeLas = [];
    for ATL = ATLs; % the different ATLs to be tested
        
        zdel = round(0.26722*fs); % delay as a fraction of the sampling rate
        
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
        n_switch = 1;
        c1_on = 1;
        
        while n<=zdel
            if n>ATL
            zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - A*zm1(n));
            zm2(n+1) = zm2(n) + f2*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + F*zm1(n-ATL));
            n = n + 1;
            else
            zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - A*zm1(n));
            zm2(n+1) = zm2(n) + f2*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2));
            n = n + 1;                
            end
        end
        
        while n < ntime
            
            if (real(zm2(n) - zm2(n-1))) < 0 && (real(zm2(n-1) - zm2(n-2))) > 0 && n>(n_switch+400) && c1_on == 1 % true everytime a cycle ends
                c1_on = 0;
                n_switch = n;                  
            end
            if (real(zm1(n) - zm1(n-1))) < 0 && (real(zm1(n-1) - zm1(n-2))) > 0 && n>(n_switch+400) && c1_on == 0 % true everytime a cycle ends
                c1_on = 1;
                n_switch = n;
            end            
            if c1_on
                c1on(n) = 1;
            end
            if ~c1_on
                c2on(n) = 1;
            end
                                    
            zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) - A*c1on(n)*zm1(n) + F*c2on(n-ATL)*zm2(n-ATL) - D*(1/f1).*(zm1(n - zdel)));           
            zm2(n+1) = zm2(n) + f2*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) - A*c2on(n)*zm2(n) + F*c1on(n-ATL)*zm1(n-ATL) - D*(1/f2).*(zm2(n - zdel)));            
            
            n = n + 1;
        end
        
        [pks, rawbeats1] = findpeaks(real(zm1));
        ATL
        plot(real(zm1))
        hold on
        plot(rawbeats1,pks,'*')
        plot(c1on)
        pause
        
        [pks, rawbeats2] = findpeaks(real(zm2));
        
        plot(real(zm2))
        hold on
        plot(rawbeats2,pks,'*')
        plot(c2on)
        pause
        
        
        
%         ROIs = 1:length(rawbeats1)-1;
%         
%         beats1 = rawbeats1(ROIs);
%         beats2 = rawbeats2(ROIs);
%         
%         LeLa = [];
%         for i=1:2:length(beats1)-1
%             LeLa = [LeLa (beats2(i) - beats1(i))];
%         end
%         for i=2:2:length(beats1)-1
%             LeLa = [LeLa (beats1(i) - beats2(i))];
%         end
%         
%         LeLas = [LeLas mean(LeLa)];
    end
%     LeLa_all = [LeLa_all; LeLas];
end

subplot(2,2,4)
errorbar(ATLs/10,100*mean(LeLa_all,1)/6660,100*std(LeLa_all,1)/6660,'o-k','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1.25)
hold on
LeLa_fit = polyfit(ATLs,100*mean(LeLa_all,1)/6660,1);
LeLa_fit = polyval(LeLa_fit,ATLs);
plot(ATLs/10,LeLa_fit,'color',[0.5 0.5 0.5],'LineWidth',2)
tx = text(0.02,0.98,'D','Units', 'Normalized', 'VerticalAlignment', 'Top');
tx.FontSize = 20;
idem = idem + 1;
ylim([-13 4])
xlim([-2 82])
grid on
set(gca,'FontSize',fontsize)
xlabel('Transmission Latency (ms)')