clear all
close all
clc

fontsize = 15;

% time
fs = 10000;
T = 1/fs;
dur = 200.1;
t = 0:T:(dur-T);
ntime = length(t);
halfsamps = floor(ntime/2);

% fixed parameters
F = 1;
f = 1;
a = 1;
b = -1;
D = 1.0;
for A = [-1.0 -0.5 0.0 0.5 1.0]
    
    
    % parameters to iterate over
    taus = linspace(0,0.5,51);
    phases = pi/2;%linspace(pi,4);
    
    for j=1:length(phases)
        biff = zeros(size(taus));
        
        for i=1:length(taus)
            
            x = exp(1i*2*pi*t*f);
            zdel= round(fs*taus(i));
            zm = ones(size(t))*exp(1i*2*pi);
            
            % main loop
            for n = 1:ntime-1
                
                if n > zdel
                    
                    zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)) - D*(1/(f)).*(zm(n - zdel)));
                    
                else
                    zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)));
                end
                
            end
            
            try
                % Peaks for oscillator and stimilus
                [pks_F,locs_F] = findpeaks(real(x));
                %             plot(real(x))
                %             hold on
                [pks_z,locs_z] = findpeaks(real(zm));
                %             plot(real(zm))
                %             plot(locs_F,pks_F,'*')
                locs_F = [1 locs_F];
                %             plot(locs_z,pks_z,'*')
                %             pause
                %             clf
                
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
                
            catch
                biff(i) = 0;
            end
            
            
            phase_diff = unwrap(angle(x)-angle(zm));
            biff(i) = mean(phase_diff((length(phase_diff)/2):end));
            
            %             plot(t,unwrap(phase_diff))
            %             xlim([0 100])
            %             ylim([-pi pi])
            %             taus(i)
            %             biff(i)
            %             pause
            %             clf
            
        end
        
        figure('position', [0, 0, 700, 300])
        plot(taus,biff,'k','LineWidth',2)
        grid on
        ylim([-pi pi])
        xlabel('\bf{\tau}','FontSize',fontsize)
        ylabel('Asynchrony (radians)','FontSize',fontsize)
        set(gca,'FontSize',fontsize)
        pause
        clf
    end
    
end

%%

clear all
close all
clc

fontsize = 15;

% time
fs = 10000;
T = 1/fs;
dur = 200.1;
t = 0:T:(dur-T);
ntime = length(t);
halfsamps = floor(ntime/2);

% fixed parameters
F = 1;
f = 1;
a = 1;
b = -1;
D = 1.0;
for A = [1.0]
    
    
    % parameters to iterate over
    taus = [0.4];
    phases = linspace(-pi,pi,9);
    phases = phases(1:end-1);
    rhos = linspace(0.25, 1.25, 3);
    
    for i=1:length(taus)
        zdel= round(fs*taus(i));
        figure('position', [0, 0, 200, 200])
        for k=1:length(rhos)
            rho = rhos(k);
            for j=1:length(phases)
                phase = phases(j);
                biff = zeros(size(taus));
                
                x = exp(1i*2*pi*t*f);
                
                zm = rho*ones(size(t))*exp(1i*phase);
                
                % main loop
                for n = 1:ntime-1
                    
                    if n > zdel
                        
                        zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)) - D*(1/(f)).*(zm(n - zdel)));
                        
                    else
                        zm(n+1) = zm(n) + f.*T*(zm(n)*(a + 1i*2*pi + b*abs(zm(n)).^2) + (F*x(n) + A*zm(n))/abs(F*x(n) + A*zm(n)));
                    end
                    
                end
                
                phase_diff = unwrap(angle(x)-angle(zm));
                amp_diff = abs(zm);
                                
%                 polarplot(phase_diff,amp_diff,'b')
                plot(real(zm))                
                hold on
                plot(real(x))
%                 A
%                 taus(i)
%                 rho
%                 phase
                pause
            end
        end
        pause
        clf
    end
    
end