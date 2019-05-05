clear all
close all
clc

fontsize = 12.5;

subplot(2,3,1) % behavioral data (solo)
NMAs = -[5.7 15; 0 0];
NMAstds = [2.4 4; 0 0];
hold on
grid on
errorbar(1.75,-5.7,2.4,'k','Marker','o','LineStyle','none','MarkerFaceColor','w','MarkerSize',10)
hold on
errorbar(2.25,-15,4,'k','Marker','o','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
Labels = {'Musician (Solo)'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
legend('Feedback On','Feedback Off','Location','southwest')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

subplot(2,3,4) % behavioral data (duet)
NMAs = -[6.8 19 ; 0 0];
NMAstds = [2 3; 0 0];
hold on
grid on
errorbar(1.75,-6.8,2,'k','Marker','o','LineStyle','none','MarkerFaceColor','w','MarkerSize',10)
hold on
errorbar(2.25,-19,3,'k','Marker','o','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
Labels = {'Musician (Duet)'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
legend('Feedback On','Feedback Off','Location','southwest')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

% musician simulations

% time paramters
fs = 10000;
T = 1/fs;
dur = 20.1;
t = 0:T:(dur-T);
ntime = length(t);

% metronome parameters
IOI = 1.0;
ft = 1/IOI;
bps = ft;
bpm = bps*60;
x = exp(1i*2*pi*t*ft);
F = 1;

% oscillator parameters
lag1cs = [];
all_NMAs = [];
a = 1;
b = -1;
D = 0.05;

zdel = round(0.222*fs); % max delay as a fraction of the sampling rate

%%% solo %%%
% oscillator initialization
zm1 = zeros(size(t));
f1 = 1;

A = -0.5;

n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        %%% player 1 plays%%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n) - D*(1/f1).*(zm1(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
locs = [1 locs];
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs;

NMA1 = zm1_locs - x_locs;
NMA = mean(NMA1((floor(length(NMA1)/2)):end)/fs)*1000;

subplot(2,3,2)
hold on
grid on
plot(1.75,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor','w','MarkerSize',10)
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

A = 0;

n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        %%% player 1 plays%%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n) - D*(1/f1).*(zm1(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
locs = [1 locs];
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs;

NMA1 = zm1_locs - x_locs;
NMA = mean(NMA1((floor(length(NMA1)/2)):end)/fs)*1000;

subplot(2,3,2)
hold on
grid on
plot(2.25,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
Labels = {'Musician SAPPA (Solo)'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

% nonmusician simulations

D = 0.36;

A = -0.5;

n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        %%% player 1 plays%%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n) - D*(1/f1).*(zm1(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
locs = [1 locs];
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs;

NMA1 = zm1_locs - x_locs;
NMA = mean(NMA1((floor(length(NMA1)/2)):end)/fs)*1000;

subplot(2,3,3)
hold on
grid on
plot(1.75,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor','w','MarkerSize',10)
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

A = 0;

n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        %%% player 1 plays%%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + F*x(n) - D*(1/f1).*(zm1(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
locs = [1 locs];
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs;

NMA1 = zm1_locs - x_locs;
NMA = mean(NMA1((floor(length(NMA1)/2)):end)/fs)*1000;

subplot(2,3,3)
hold on
grid on
plot(2.25,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
Labels = {'Non-musician SAPPA (Solo)'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

%%% duet simulations %%%

% musicians

D = 0.05;

A = -0.5;

zm1 = zeros(size(t));
zm2 = zeros(size(t));

n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) - A*zm2(n)));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) - A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm2(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm2(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs(2:2:end);

NMA1 = zm1_locs - x_locs;

% find the time of metronome beats
[pks,locs] = findpeaks(-real(x));
x_locs = locs(2:end);

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm2));
zm1_locs = locs(2:2:end);

NMA2 = zm1_locs - x_locs;

NMA = mean([NMA1((floor(length(NMA1)/2)):end) NMA2((floor(length(NMA2)/2)):end)]/fs)*1000;

subplot(2,3,5)
hold on
grid on
plot(1.75,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor','w','MarkerSize',10)
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

A = 0;
n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) - A*zm2(n)));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) - A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm2(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm2(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
locs = [1 locs];
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs;

NMA1 = zm1_locs - x_locs;

% find the time of metronome beats
[pks,locs] = findpeaks(-real(x));
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm2));
zm1_locs = locs;

NMA2 = zm1_locs - x_locs;

NMA = mean([NMA1((floor(length(NMA1)/2)):end) NMA2((floor(length(NMA2)/2)):end)]/fs)*1000;

subplot(2,3,5)
hold on
grid on
plot(2.25,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
Labels = {'Musician SAPPA (Duet)'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

% non-musicians

D = 0.36;

A = -0.5;

zm1 = zeros(size(t));
zm2 = zeros(size(t));

n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) - A*zm2(n)));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) - A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm2(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm2(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
x_locs = [1 locs];

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs;

NMA1 = zm1_locs - x_locs;

% find the time of metronome beats
[pks,locs] = findpeaks(-real(x));
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm2));
zm1_locs = locs;

NMA2 = zm1_locs - x_locs;

NMA = mean([NMA1((floor(length(NMA1)/2)):end) NMA2((floor(length(NMA2)/2)):end)]/fs)*1000;

subplot(2,3,6)
hold on
grid on
plot(1.75,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor','w','MarkerSize',10)
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')

A = 0;
n = 1;
while n<ntime
    
    if n<=zdel
        %%% delay buffer %%%
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) - A*zm2(n)));
    elseif mod(t(n)*10,10) > 5 && mod(t(n)*10,10) < 10
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) + A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) - A*zm1(n))/abs(F*x(n) + A*zm1(n)) - D*(1/f1).*(zm2(n - zdel)));
    else
        zm1(n+1) = zm1(n) + f1*T*(zm1(n)*(a + 1i*2*pi + b*abs(zm1(n)).^2) + (F*x(n) - A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm1(n - zdel)));
        zm2(n+1) = zm2(n) + f1*T*(zm2(n)*(a + 1i*2*pi + b*abs(zm2(n)).^2) + (-F*x(n) + A*zm2(n))/abs(F*x(n) + A*zm2(n)) - D*(1/f1).*(zm2(n - zdel)));
    end
    n = n + 1;
end

% find the time of metronome beats
[pks,locs] = findpeaks(real(x));
locs = [1 locs];
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm1));
zm1_locs = locs;

NMA1 = zm1_locs - x_locs;

% find the time of metronome beats
[pks,locs] = findpeaks(-real(x));
x_locs = locs;

% find the time of zm1 beats
[pks,locs] = findpeaks(real(zm2));
zm1_locs = locs;

NMA2 = zm1_locs - x_locs;

NMA = mean([NMA1((floor(length(NMA1)/2)):end) NMA2((floor(length(NMA2)/2)):end)]/fs)*1000;

subplot(2,3,6)
hold on
grid on
plot(2.25,NMA,'k','Marker','o','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
Labels = {'Non-musician SAPPA (Duet)'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylim([-80.5 0])
xlim([1.2 2.8])
ylabel('Mean Asynchrony (ms)')
set(gca,'FontSize',fontsize)
set(gca,'xaxisLocation','top')