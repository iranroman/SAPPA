% clear all
close all
clc

fontsize = 15;

% time paramters
fs = 10000;
T = 1/fs;
dur = 100.1; % seconds
t = 0:T:(dur+T);
ntime = length(t);
halfsamps = floor(ntime/2);

% stimulus
ft = 1.0;
x = ones(size(t));
x(mod(1:ntime,fs*(1/ft))>fs*(1/ft)/2) = 0;
x = 2*x-1;
        
% oscillator initialization and parameters
a = 1;
b = -1;
zm = 0.5*ones(size(t))*exp(1i*2*pi); % SMS oscillation
f = ft;

% parameters we will iterate over
As = [-0.5];
Fs = 1;
Ds = linspace(0,0.4,10);
taus = linspace(0,0.3,10);

for iA = 1:length(As)
    
    A = As(iA);
    
    F = Fs;
    
    X = NaN(length(Ds),length(taus));
    
    for iD = 1:length(Ds)
        
        D = Ds(iD);
        
        for itau = 1:length(taus)
            
            tau = taus(itau);
            
            zdel = round(tau*fs); % SMS delay as a fraction of the sampling rate
            
            % forward euler integration loop
            synch = 0;
            for n = 1:ntime-1
                
                    if n > zdel && synch == 0
                        zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)) - D*(1/(f)).*(zm(n - zdel)));
                    else
                        zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)));
                    end

            end
            
            if synch == 0
                % Peaks for oscillator and stimilus
                [pks_F,locs_F] = findpeaks(sin(2*pi*t*ft));
                %             plot(real(x))
                %             hold on
                [pks_z,locs_z] = findpeaks(real(zm));
                %             plot(real(zm))
                %             plot(locs_F,pks_F,'*')
                locs_F = [1 locs_F];
                %             plot(locs_z,pks_z,'*')
                %             pause
                %             clf
                
                
                
                try
                    % which z peak is closest to the midpoint of the
                    % simulation?
                    halfsamps_locsz_diff = abs(halfsamps-locs_z);
                    [~,mid_nzpeak_index] = min(halfsamps_locsz_diff);
                    mid_nzpeak = locs_z(mid_nzpeak_index);
                    
                    % eliminate the first half of the simulation for z
                    locs_z = locs_z(mid_nzpeak_index:end);
                    
                    % which F peak is closest to mid_nzpeak?
                    mid_nzpeak_locs_F_diff = abs(locs_F - mid_nzpeak);
                    [~,mid_F_peaks_index] = min(mid_nzpeak_locs_F_diff);
                    
                    % which z peak is the penultimate one?
                    pen_nzpeak = locs_z(end-1);
                    % which F peak is closest to the penultimate z peak?
                    pen_nzpeak_locs_F_diff = abs(locs_F - pen_nzpeak);
                    [~,pen_F_peaks_index] = min(pen_nzpeak_locs_F_diff);
                    
                    % compute the mean asynchrony
                    mean_asynchrony = locs_z(1:end-1) - locs_F(mid_F_peaks_index:pen_F_peaks_index);                                                            
                    
                    X(iD,itau) = 1000*mean(mean_asynchrony)/fs;
                    
                catch
                    
                end
            else
                X(iD,itau) = 0;
            end
%             X(iD,itau)
%             plot(real(x))
%             hold on
%             plot(real(zm))
%             pause
%             clf
        end
    end
    
%     subplot(1,3,2)
    imagesc(X)
    xlabels = zeros(1,length(taus)*2+1);
    xlabels(1,2:2:end) = taus;
    xlabels = num2cell(xlabels);
    xlabels(1,1:2:end) = {['']};
    ylabels = zeros(1,length(Ds)*2+1);
    ylabels(1,2:2:end) = Ds;
    ylabels = num2cell(ylabels);
    ylabels(1,1:2:end) = {['']};
    xticklabels(round(taus,2))
    yticklabels(round(Ds,2))
    title(sprintf('A=%.2f',A))
    xlabel('\bf{\tau}')
    ylabel('D')
    caxis([-300 10])
%     c = colorbar;
%     c.Label.String = 'Asynchrony';
    colormap(parula)
    set(gca,'FontSize',fontsize)
    if true
        % code
        [m n]=size(X);
        hold on;
        for i = 1:m
            for j = 1:n
                nu = X(i,j);
                val = num2str(round(nu));
                text(j,i,val,'FontSize',13)
            end
        end
        hold off;
    end
        
end
%%
%%% a sketch explaining how the NMA is computed %%%
fs = 1000;
T = 1/fs;
dur = 0.85;
t = -0.85:T:(dur+T);
ntime = length(t);

x = ones(size(t));
x(mod(401:ntime+400,fs*(1/ft))>fs*(1/ft)/2) = 0;
x = 2*x-1;
zm = 0.6.*cos(2*pi*t+pi/6);

% find the peaks
[pks, rawbeats1] = findpeaks(real(x));
[pks, rawbeats2] = findpeaks(real(zm));
rawbeats1 = rawbeats1+fs/4;

subplot(1,3,1)
plot(t,x,'LineWidth',2,'DisplayName','stimulus')
hold on
plot(t,zm,'LineWidth',2,'DisplayName','agent z_m')
legend('show','Location','southeast')
plot([t(rawbeats1)-T t(rawbeats1) t(rawbeats1)+T],[-10 1 10],'--k','LineWidth',1.5)
plot([t(rawbeats2)-T t(rawbeats2) t(rawbeats2)+T],[-10 1 10],'--','Color',[0.5 0.5 0.5],'LineWidth',1)
quiver(t(rawbeats2), 1.5, 0, max(real(zm)) + 0.1 - 1.5,'k','LineWidth',3,'MaxHeadSize',1)
tx = text(0.2,0.85,'Asynchrony','Units', 'Normalized', 'VerticalAlignment', 'Top','BackgroundColor','white','EdgeColor','k','LineWidth',1);
tx.FontSize = fontsize;
xlim([min(t) max(t)])
ylim([-2.5 2.5])
ylabel('Amplitude')
xlabel('time (s)')
set(gca,'FontSize',fontsize)
grid on
%%
fontsize = 15;
% time paramters
fs = 1000;
T = 1/fs;
dur = 200; % seconds
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

for D = [0.04 0.14];
    NMAs = [];
    NMAstds = [];
    for IOI = IOIs; % in seconds
        ft = 1/IOI;
        bps = ft;
        bpm = bps*60;
        x = cos(2*pi*t*ft);
        x = ones(size(t));
        x(mod(1:ntime,fs*(1/ft))>fs*(1/ft)/2) = 0;
        x = 2*x-1;
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
        
%         IOI
%         plot(t,sin(2*pi*t*ft))
%         hold on
%         plot(t,real(zm))
%         ylim([-1.5 1.5])
%         pause
%         clf
        % finding the peaks and the NMA
        [pks, rawbeats1] = findpeaks(sin(2*pi*t*ft));
        [pks, rawbeats2] = findpeaks(real(zm));
%         rawbeats1 = [1 rawbeats1]; % adding the peak at n  = 1
        NMA = rawbeats2 - rawbeats1;
        NMAs = [NMAs mean(NMA(floor(end/2):end))];
        
    end
    
    subplot(1,3,3)
    plot(1000*IOIs, 1000*NMAs/fs,'o','MarkerSize',10,'LineWidth',2,'DisplayName',sprintf('D=%.3f (A=-0.5)',D),'color',cols{icolor},'MarkerFaceColor',cols{icolor})        
    grid on
    grid minor
    hold on
    xlim([1000*IOIs(1) - 100 1000*IOIs(end) + 100])
    ylim([-170 20])
    xlabel('Stimulus Period (ms)','FontSize',fontsize)
    ylabel('Asynchrony (ms)','FontSize',fontsize)
    set(gca,'FontSize',fontsize)
    icolor = icolor + 1;
end

nmNMA = 1000*[-0.036 -0.048 -0.051 -0.066 -0.0775 -0.077 -0.103 -0.111 -0.1065 -0.142 -0.1475]; % From Repp 2007
nmNMAse = 1000*[0.008 0.012 0.013 0.015 0.016 0.015 0.013 0.021 0.022 0.025 0.030];
mNMA = 1000*[-0.0045 -0.008 -0.014 -0.014 -0.003 -0.007 0.000 -0.009 -0.006 -0.031 -0.031];
mNMAse = 1000*[0.0025 0.003 0.0055 0.008 0.008 0.007 0.011 0.0085 0.016 0.0125 0.016];

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