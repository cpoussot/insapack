clc, clearvars, close all, 
%
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex');
end
% MDSPACK
setenv('MDSHOME','/Users/charles/Documents/MDS')
addpath('/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/v1.1.0/API/matlab/')
addpath('/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/v1.1.0/bin')
%
addpath('/Users/charles/Mon Drive/Research/Codes/')
%
col     = colororder;
fmax    = 20;
fBnd    = [0 fmax];
Ts      = 1/fmax/5;
Ns      = 2^10+1;
FLIP    = false;
RPHI    = true;
P       = 20;
noise_i = 1e-3;
noise_o = 2e-2;
%
rng(1)
%G       = tf(1,[1/(2*pi)^2 2*.01/2*pi 1]);
%K       = pidtune(G,'pi',1);
%
G       = stabsep(rss(200,1,1)); 
K       = pidtune(G,'pi',.1);
CL      = [feedback(G*K,1); K/(1+K*G)];

kk      = 0;
opt     = struct('f_band',fBnd,'flip',FLIP,'r_phi',RPHI);
[r,t]   = mdspack.signals('MSIN',Ns,Ts,opt);
%[r,t]   = mdspack.signals('MLBS',Ns,Ts,opt);
%[r,t]   = mdspack.signals('CHRP',Ns,Ts,opt);

out = lsim(CL,r,t);
y   = out(:,1);
u   = out(:,2);

% for i = 1:P
%     % noise
%     nin     = noise_i*randn(nt,1);
%     ncon    = noise_i*randn(nt,1);
%     nout    = noise_o*randn(nt,1);
%     % i/o
%     rn(:,i) = r.*(1+nin);
%     un(:,i) = u.*(1+ncon);
%     yn(:,i) = y.*(1+nout)+0.*t;
% end

nt  = length(t);
f   = linspace(0,1,nt).'/Ts;
R   = fft(r)/nt;
U   = fft(u)/nt;
Y   = fft(y)/nt;

SYR = sum((Y.*conj(R)),2)./sum((R.^2),2);
SUR = sum((U.*conj(R)),2)./sum((R.^2),2);

Ge  = SYR./SUR;
Gt  = freqresp(G,2*pi*f);
Gcl = freqresp(CL(1,1),2*pi*f);
%
figure
%
kk = kk+1; subplot(3,2,kk); hold on, axis tight,
title('Reference signal'), xlabel('Time (s)'), ylabel('Amplitude')
plot(t,r,'k-'), 
%
kk = kk+1; subplot(3,2,kk); hold on, axis tight, 
plot(t,u,'k-')
title('Control signal'), xlabel('Time (s)'), ylabel('Amplitude')
%
kk = kk+1; subplot(3,2,kk); hold on, axis tight, 
plot(t,y,'k-')
title('Output signal'), xlabel('Time (s)'), ylabel('Amplitude')
%
kk = kk+1; subplot(3,2,kk); hold on, axis tight,
title('TF gain'), xlabel('Frequency (Hz)'), ylabel('Gain (dB)')
plot(f,20*log10(abs(Gt(:))),'k-')
plot(f,20*log10(abs(Ge(:))),'r--')
plot(f,20*log10(abs(Gcl(:))),':')
set(gca,'xlim',[0 fmax],'XScale','log')
%
kk = kk+1; subplot(3,2,kk); hold on, axis tight,
title('TF phase'), xlabel('Frequency (Hz)'), ylabel('Phase (rad)')
plot(f,angle(Gt(:)),'k-')
plot(f,angle(Ge(:)),'r--')
plot(f,angle(Gcl(:)),':')
set(gca,'xlim',[0 fmax],'XScale','log')
legend({'True','Estimated (shouck inside)'})
% %
% sgtitle({['Noise i/o [\%]: $' num2str(noise_i) '~\&~ ' num2str(noise_o) '$']; ...
%          ['Periods: $' num2str(P) '$']},'FontSize',20)

