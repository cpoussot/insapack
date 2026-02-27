clear all, close all, clc;
set(groot,'DefaultFigurePosition', [100 100 1300 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for iloe1 = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(iloe1)},'factory','default'),'latex'); end

%%% INSAPACK
addpath('/Users/charles/Documents/GIT/insapack')
col     = colororder;

%%% System to identify
CAS = 3;
switch CAS
    case 1
        G       = tf([1 -2],[.1 .4 1]); G = G/dcgain(G);
        % Frequency
        w       = logspace(-2,2.5,300);
        ratio   = 5;
        maxEig  = max(abs(eig(G)));
        wmax    = ratio*maxEig;
        fmax    = wmax/2/pi;
        Fs      = 2^(nextpow2(fmax)+1);
        Ws      = 2*pi*Fs;
        Ts      = 1/Fs;
        % Duration
        Ns      = 2^9;
        % Non-parametric
        P       = 50;
        noise_i = .1;
        noise_o = .2;
        % Identification
        nx      = 2;
    case 2
        rng(1); 
        G       = stabsep(rss(20,1,1));
        % Frequency
        w       = logspace(-2,2.5,300);
        ratio   = 5;
        maxEig  = max(abs(eig(G)));
        wmax    = ratio*maxEig;
        fmax    = wmax/2/pi;
        Fs      = 2^(nextpow2(fmax)+1);
        Ws      = 2*pi*Fs;
        Ts      = 1/Fs;
        % Duration
        Ns      = 2^12;
        % Non-parametric
        P       = 50;
        noise_i = .1;
        noise_o = .2;
        % Identification
        nx      = 8%10;
    case 3
        rng(1); 
        G       = stabsep(rss(100,1,1));
        % Frequency
        w       = logspace(-2,2,300);
        ratio   = 5;
        maxEig  = max(abs(eig(G)));
        wmax    = ratio*maxEig;
        fmax    = wmax/2/pi;
        Fs      = 2^(nextpow2(fmax)+1);
        Ws      = 2*pi*Fs;
        Ts      = 1/Fs;
        % Duration
        Ns      = 2^11;
        % Non-parametric
        P       = 7;
        noise_i = .1;
        noise_o = .3;
        % Identification
        nx      = 5;
end
% %%% Noise generator
% Gn      = tf(n*[1/20 1],[1/100 1]);
% eigGn   = eig(Gn);
% frGn    = freqresp(Gn,w);
%%% Original system analysis
eigG    = eig(G);
frG     = freqresp(G,w);
theta   = linspace(pi/2,3*pi/2,500);
figure
subplot(2,2,[1 3]), 
plot(real(eigG),imag(eigG),'.','DisplayName','$\lambda(\mathbf G)$'), grid on, xlabel('Real'), ylabel('Imag.'), hold on
%plot(real(eigGn),imag(eigGn),'x','DisplayName','$\lambda(\mathbf G_n)$'), 
plot(maxEig*cos(theta),maxEig*sin(theta),'--','DisplayName','$|\lambda_{max}|$'),
plot(wmax*cos(theta),wmax*sin(theta),'--','DisplayName',['$' num2str(ratio) '|\lambda_{max}|$']),
plot(Ws*cos(theta),Ws*sin(theta),'--','DisplayName','$2\pi F_s$'),
hh = gca;
plot([0 0],hh.YLim,'k:','DisplayName','Stability limit')
legend('show','location','best'), axis equal
title('Eigenvalues')
subplot(222)
plot(w,20*log10(abs(frG(:))),'-','DisplayName','$\mathbf G$'), set(gca,'XScale','log'), grid on, xlabel('Pulsation [rad/s]'), ylabel('Gain [dB]'), hold on
%plot(w,20*log10(abs(frGn(:))),'-','DisplayName','$\mathbf G_n$'), set(gca,'XScale','log'), grid on, xlabel('Pulsation [rad/s]'), ylabel('Gain [dB]'), hold on
axis tight, hh = gca;
plot([1 1]*maxEig,hh.YLim,'--','DisplayName','$|\lambda_{max}|$'),
plot([1 1]*wmax,hh.YLim,'--','DisplayName',['$' num2str(ratio) '|\lambda_{max}|$']),
plot([1 1]*Ws,hh.YLim,'--','DisplayName','$2\pi F_s$')
legend('show','location','best')
title('Bode gain')
subplot(224)
plot(w,angle(frG(:)),'-','DisplayName','$\mathbf G$'), set(gca,'XScale','log'), grid on, xlabel('Pulsation [rad/s]'), ylabel('Gain [dB]'), hold on
%plot(w,angl(frGn(:)),'-','DisplayName','$\mathbf G_n$'), set(gca,'XScale','log'), grid on, xlabel('Pulsation [rad/s]'), ylabel('Gain [dB]'), hold on
axis tight, hh = gca;
plot([1 1]*maxEig,hh.YLim,'--','DisplayName','$|\lambda_{max}|$'),
plot([1 1]*wmax,hh.YLim,'--','DisplayName',['$' num2str(ratio) '|\lambda_{max}|$']),
plot([1 1]*Ws,hh.YLim,'--','DisplayName','$2\pi F_s$')
legend('show','location','best')
title('Bode phase')
%%
%%% Generate exciting signal
FBND        = [0 Fs/4];
REV         = false;
SHOW        = true;
RPHI        = false;
ODD         = 'all';
[u,t,info]  = insapack.multisine(Ns,Ts,FBND,RPHI,ODD,REV,SHOW);
u = u.'; t = t.';

%%% Apply to system to be identified
y   = lsim(G,u,t);
for i = 1:P
    % +i/o noise
    nin     = noise_i*randn(length(t),1);
    nout    = noise_o*randn(length(t),1);
    un(:,i) = u.*(1+nin);
    yn(:,i) = y.*(1+nout);
end
[f0,U0,Y0,G0,sigU2,sigY2,sigUY2,sigG2,b,rho] = insapack.non_param_freq(un,yn,Ts);
w0 = 2*pi*f0;

%%% Data
figure, 
subplot(221); hold on, grid on, axis tight
h1=plot(t,yn,'-','Color',[1 1 1]*.8,'DisplayName','$\mathbf G+n$');
h2=plot(t,y,'-','Color',col(1,:),'DisplayName','$\mathbf G$');
legend('show',[h1(1) h2(1)]), 
title('Data vs. model') 
xlabel('time [s]'), ylabel('Output'),
%
subplot(222); hold on, grid on, axis tight
plot(w,20*log10(abs(frG(:))),'-','DisplayName','$\mathbf G(\imath\omega)$')
plot(w0,20*log10(abs(G0(:))),'.','DisplayName','$\mathbf G_0(\imath\omega)$')
hh = gca;
plot([1 1]*wmax,hh.YLim,'--','DisplayName',['$' num2str(ratio) '|\lambda_{max}|$']),
plot([1 1]*Ws/2,hh.YLim,'--','DisplayName','$\pi F_s$')
plot([1 1]*2*pi*max(FBND),hh.YLim,'--','DisplayName','Max. exct. signal')
set(gca,'XScale','log')
legend('show','Location','best')
title('Bode gain')
xlabel('Pulsation [rad/s]'), ylabel('Gain [dB]'),
%
subplot(224); hold on, grid on, axis tight
plot(w,angle(frG(:)),'-','DisplayName','$\mathbf  G(\imath\omega)$')
plot(w0,angle(G0(:)),'.','DisplayName','$\mathbf G_0(\imath\omega)$')
hh = gca;
plot([1 1]*wmax,hh.YLim,'--','DisplayName',['$' num2str(ratio) '|\lambda_{max}|$']),
plot([1 1]*Ws/2,hh.YLim,'--','DisplayName','$\pi F_s$')
plot([1 1]*2*pi*max(FBND),hh.YLim,'--','DisplayName','Max. exct. signal')
set(gca,'XScale','log')
legend('show','Location','best')
title('Bode phase')
xlabel('Pulsation [rad/s]'), ylabel('Phase [rad]'),
%%
%%% Identification via N4SID
% Hn4sid          = n4sid(u,mean(yn,2),nx,'Ts',Ts);
% Hn4sid          = d2c(stabsep(ss(Hn4sid)),'tustin');
% frHn4sid        = freqresp(Hn4sid,w);
% eigHn4sid       = eig(Hn4sid);

%%% Identification via N4SID in frequency-domain
data            = iddata(Y0,U0,Ts,'Frequency',w0)
Hn4sid          = n4sid(data,nx);
Hn4sid          = d2c(stabsep(ss(Hn4sid)),'tustin');
frHn4sid        = freqresp(Hn4sid,w);
eigHn4sid       = eig(Hn4sid);

%%% Identification via Loewner
wid             = 2*pi*f0(f0<FBND(end)/2);
wRange          = 1:floor(length(wid)/2)*2;
wid             = 2*pi*f0(wRange);
G0id            = G0(wRange);
[la,mu,W,V,R,L] = insapack.data2loewner(wid,G0);
opt.target      = nx;
[hr,info]       = insapack.loewner_tng(la,mu,W,V,R,L,opt);
Hloe            = dss(info.Ar,info.Br,info.Cr,info.Dr,info.Er);
Hloe            = stabsep(Hloe);
frHloe          = freqresp(Hloe,w);
eigHloe         = eig(Hloe);

%%% Validation signal
[uv,tv,infov]   = insapack.mlbs(Ns,Ts,FBND,REV,SHOW); uv = uv.'; tv = tv.';
yv              = lsim(G,uv,tv);
yn4sid          = lsim(Hn4sid,uv,tv);
yloe            = lsim(Hloe,uv,tv);

figure
subplot(221); hold on, grid on, axis tight
plot(tv,yv,'-','DisplayName','$\mathbf G$');
plot(tv,yn4sid,'--','DisplayName','$\mathbf H_{n4sid}$');
plot(tv,yloe,'--','DisplayName','$\mathbf H_{loe}$');
legend('show'), 
subplot(223),hold on, grid on, axis equal
plot(real(eigG),imag(eigG),'o','DisplayName','$\lambda(\mathbf G)$'), grid on, 
plot(real(eigHn4sid),imag(eigHn4sid),'x','DisplayName','$\lambda(\mathbf H_{n4sid})$')
plot(real(eigHloe),imag(eigHloe),'s','DisplayName','$\lambda(\mathbf H_{loe})$')
plot(maxEig*cos(theta),maxEig*sin(theta),'k--','DisplayName','$|\lambda_{max}|$'),
hh = gca;
plot([0 0],hh.YLim,'k:','DisplayName','Stability limit')
legend('show','Location','best'), 
xlabel('Real'), ylabel('Imag.'), hold on
subplot(222); hold on, grid on, axis tight
plot(w,20*log10(abs(frG(:))),'-','DisplayName','$\mathbf G(\imath\omega)$')
plot(w,20*log10(abs(frHn4sid(:))),'--','DisplayName','$\mathbf H_{n4sid}$')
plot(w,20*log10(abs(frHloe(:))),'--','DisplayName','$\mathbf H_{loe}$')
set(gca,'XScale','log'), 
legend('show','Location','best')
title('Bode gain'), xlabel('Pulsation [rad/s]'), ylabel('Gain [dB]'),
subplot(224); hold on, grid on, axis tight
plot(w,angle(frG(:)),'-','DisplayName','$\mathbf G(\imath\omega)$')
plot(w,angle(frHn4sid(:)),'--','DisplayName','$\mathbf H_{n4sid}$')
plot(w,angle(frHloe(:)),'--','DisplayName','$\mathbf H_{loe}$')
set(gca,'XScale','log'), 
legend('show','Location','best')
title('Bode phase'), xlabel('Pulsation [rad/s]'), ylabel('Phase [rad]'),

%%% Metrics
N       = length(yv);
En4sid  = 1/N*sum(abs(yv-yn4sid)/max(abs(yv)));
Eloe    = 1/N*sum(abs(yv-yloe)/max(abs(yv)));
sgtitle(sprintf('TD errors: $%0.2f$ (N2SID) / $%0.2f$ (LF)',En4sid,Eloe),'interpreter','latex','FontSize',20) 
% %%
% figure, 
% plot(info.sv,'-o'), grid on, axis tight
% set(gca,'YScale','log')
% xlabel('Index $k$'), ylabel('Singular value'), title('Normalized Loewner singular value')