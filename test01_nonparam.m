clc, clearvars, close all, 
%
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end;
%
col     = colororder;
% Signal
fmax    = 20;
fBnd    = [0 fmax];
Ts      = 1/fmax/10;
Ns      = 2^12+1;
Res     = 1/(Ns*Ts) % frequency resolution
RPHI    = true;
ODD     = false;
FLIP    = false;
SHOW    = false;
P       = 20;
% System
noise_i = 1e-1;
noise_o = 2e-1;
biais_o = 0;
G       = stabsep(rss(200,1,1)); 
%G       = freqsep(G,fBnd*2*pi);

kk          = 0;
[u,t,info]  = insapack.multisine(Ns,Ts,fBnd,RPHI,ODD,FLIP,SHOW); u = u.'; t = t.';

y   = lsim(G,u,t);
nt  = length(t);
U   = fft(u)/nt;
Y   = fft(y)/nt;
for i = 1:P
    % noise
    nin     = noise_i*randn(nt,1);
    nout    = noise_o*randn(nt,1);
    % i/o
    un(:,i) = u.*(1+nin);
    yn(:,i) = y.*(1+nout)+biais_o*randn(1).*t;
end
[f,U0,Y0,G0,sigU2,sigY2,sigUY2,sigG2,b,rho] = insapack.non_param_freq(un,yn,Ts);
%
figure
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight,
title('Input signal'), xlabel('Time (s)'), ylabel('Amplitude')
plot(t,un,'-','Color',[1 1 1]*.8,'LineWidth',1)
plot(t,u,'k-'), 
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight, 
plot(t,yn,'-','Color',[1 1 1]*.8,'LineWidth',1)
plot(t,y,'k-')
title('Output signal'), xlabel('Time (s)'), ylabel('Amplitude')
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight,
title('Input signal'), xlabel('Frequency (Hz)'), ylabel('Amplitude')
%plot(f,20*log10(abs(Un)),'.','Color',[1 1 1]*.8)
plot(f,20*log10(abs(U)),'k-')
plot(f,20*log10(abs(U0)),'r--')
plot(f,20*log10(sigU2),'-')
set(gca,'xlim',[0 1/Ts/2],'Xscale','log')
%legend({'$|U|$ (true)','$|U_0|$ (mean)','$\sigma_U^2$'},'location','best')
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight,
title('Output signal'), xlabel('Frequency (Hz)'), ylabel('Amplitude')
%plot(f,20*log10(abs(ynf)),'.')
plot(f,20*log10(abs(Y)),'k-')
plot(f,20*log10(abs(Y0)),'r--')
plot(f,20*log10(sigY2),'-')
set(gca,'xlim',[0 1/Ts/2],'XScale','log')
%legend({'$|Y|$ (true)','$|Y_0|$ (mean)','$\sigma_Y^2$'},'location','best')
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight,
title('i/o noise correlation'), xlabel('Frequency (Hz)'), ylabel('Amplitude')
plot(f,abs(rho),'Color',col(1,:))
set(gca,'xlim',[0 1/Ts/2],'XScale','log','YLim',[0 1])
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight,
plot([0 0],[0 0],'k-')
plot([0 0],[0 0],'r--')
plot([0 0],[0 0],'-')
plot([0 0],[0 0],':','Color',col(1,:))
plot([0 0],[0 0],':','Color',[1 1 1]*.8,'LineWidth',1)
plot([0 0],[0 0],'-','Color',col(5,:))
set(gca,'XTick',[],'YTick',[])
legend({'true','Schoukens','$\sigma_U^2,\sigma_Y^2,\sigma_G^2$','\texttt{mdspack.tfest}','raw data coll.','\texttt{mdspack.loewner}'},'location','best')
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight,
% >> True
FR          = freqresp(G,2*pi*f);
% >> Estimated
[FRe,ff,~]  = tfestimate(mean(un,2),mean(yn,2),1/Ts);
% >> Identified


%
title('TF gain'), xlabel('Frequency (Hz)'), ylabel('Gain (dB)')
plot(f,20*log10(abs(FR(:))),'k-')
plot(ff,20*log10(abs(FRe(:))),'r--')
plot(f,20*log10(abs(sigG2(:))),'-')
plot(f,20*log10(abs(G0)),':','Color',col(1,:))
set(gca,'xlim',[0 1.0*fmax],'XScale','log')
%
kk = kk+1; subplot(4,2,kk); hold on, axis tight,
title('TF phase'), xlabel('Frequency (Hz)'), ylabel('Phase (rad)')
plot(f,angle(FR(:)),'k-')
plot(ff,angle(FRe(:)),'r--')
plot(f,angle(Y0./U0),':','Color',col(1,:))
set(gca,'xlim',[0 1.0*fmax],'XScale','log')
%
sgtitle({['Noise i/o [\%]: $' num2str(100*noise_i) '~\&~ ' num2str(100*noise_o) '$ and biais ' num2str(biais_o)]; ...
         ['Periods: $' num2str(P) '$']},'FontSize',20)

