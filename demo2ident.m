% DEMO IDENTIFICATION 
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
% Date  : March 2022 (creation)
%         February 2026 (modification)
% 
% Description
% Demonstration script for data-driven dynamical model identification using 
% either
% 
% Note
% The script uses the "+insa" matlab package.
%
clearvars, close all, clc
set(groot,'DefaultFigurePosition', [300 100 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for iloe1 = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(iloe1)},'factory','default'),'latex'); end
clear all, close all, clc;

%%% INSAPACK & LF
addpath('/Users/charles/Documents/GIT/insapack')
addpath('/Users/charles/Documents/GIT/lf')

%%% LOAD DATA
load('signals')
sig     = sig2;
Ts      = sig.Ts;
t       = sig.t;
u       = sig.u;
y       = sig.y;
yn      = sig.yn;
ytrue   = lsim(G,u,t);
w       = logspace(-2,log10(1/Ts),200)*2*pi;
Gtrue   = freqresp(G,w);

%%% IDENTIFICATION VIA IDENT
nx      = 3;
Hid     = n4sid(u,yn,nx,'Ts',Ts);
%
yid     = lsim(Hid,u,t);
Gid     = freqresp(Hid,w);

%%% IDENTIFICATION VIA LOEWNER
[f,U0,Y0,G0]    = insapack.non_param_freq(u,yn,Ts);
puls            = 2*pi*f(f<FBND(end)/2);
wRange          = 1:floor(length(puls)/2)*2;
puls            = 2*pi*f(wRange);
G0              = G0(wRange);
%
[la,mu,W,V,R,L] = insapack.data2loewner(puls,G0);
opt             = [];
opt.target      = nx;
[hr,info]       = lf.loewner_tng(la,mu,W,V,R,L,opt);
%
Hloe            = dss(info.Ar,info.Br,info.Cr,info.Dr,info.Er);
Hloe            = stabsep(Hloe);
yloe            = lsim(Hloe,u,t);
Gloe            = freqresp(Hloe,w);

%%
%%% CHECK
figure, 
subplot(221); hold on, grid on, axis tight
plot(t,ytrue,'-','DisplayName','$G$')
plot(t,yid,'--','DisplayName','$H_{n4sid}$')
plot(t,yloe,'--','DisplayName','$H_{loe}$')
%plot(t,yn,':','DisplayName','$y_G(kT_s)+n$')
legend('show'), set(gca,'YLim',[min(ytrue) max(ytrue)])
title('Data vs. model')
%
subplot(222); hold on, grid on, axis tight
plot(w,20*log10(abs(Gtrue(:))),'-','DisplayName','$G(\imath\omega)$')
plot(w,20*log10(abs(Gid(:))),'--','DisplayName','$H_{n4sid}(\imath\omega)$')
plot(w,20*log10(abs(Gloe(:))),'--','DisplayName','$H_{loe}(\imath\omega)$')
plot(puls,20*log10(abs(G0(:))),':','DisplayName','$\hat G(\imath\omega)$')
set(gca,'XScale','log')
legend('show','Location','best')
title(['Bode gain (Nyquist, $' num2str(2*pi*1/Ts) '$)'])
subplot(224); hold on, grid on, axis tight
plot(w,angle(Gtrue(:)),'-','DisplayName','$G(\imath\omega)$')
plot(w,angle(Gid(:)),'--','DisplayName','$H_{n4sid}(\imath\omega)$')
plot(w,angle(Gloe(:)),'--','DisplayName','$H_{loe}(\imath\omega)$')
plot(puls,angle(G0(:)),':','DisplayName','$\hat G(\imath\omega)$')
set(gca,'XScale','log')
legend('show','Location','best')
title(['Bode phase (Nyquist, $' num2str(2*pi*1/Ts) '$)'])

