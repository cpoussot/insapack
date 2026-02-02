% DEMO IDENTIFICATION 
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
% Date  : March 2022 (creation)
% 
% Description
% Demonstration script for data-driven dynamical model identification using 
% either
%  * the "ident" toolbox embedded in the Matlab Identification Toolbox, or 
%  * the Hankel driven 
% This script is a rather simple one and modifications may be applied by 
% users. Still, it provides a first shot for it, but does not prevent for 
% additional reading and thinking.
% 
% Note
% The script uses the "+insa" matlab package.
% 
% setenv('MDSHOME','/Users/charles/Documents/MDS')
% addpath('/Users/charles/Documents/MDS/mdspack/API/matlab/')
% addpath('/Users/charles/Documents/MDS/mdspack/build/OSX/lib')

clear all, close all, clc;
% >> Graphic settings
insa.initGraphics
load('signals')

% Identify using 'ident' on either (u,y) or (u,yn) couples
t1  = t1(:); u1  = u1(:); y1  = y1n(:); dt1 = t1(2)-t1(1);
t2  = t2(:); u2  = u2(:); y2  = y2n(:); dt2 = t2(2)-t2(1);
t3  = t3(:); u3  = u3(:); y3  = y3n(:); dt3 = t3(2)-t3(1);
t4  = t4(:); u4  = u4(:); y4  = y4n(:); dt4 = t4(2)-t4(1);

%%% >> IDENTIFICATION VIA IDENT
systemIdentification

% %%% >> IDENTIFICATION VIA HANKEL (JUST FOR FUN)
% load('signals')
% tt = t3(:); uu = u3(:); yy = y3(:); dt = tt(2)-tt(1);
% r       = [];
% Hr      = insa.pencil(uu,yy,dt,r);
% Hr      = stabsep(Hr);
% Hr      = minreal(Hr);
% Hr      = d2c(Hr,'zoh');
% size(Hr)
% % Time-domain simulation
% tt2     = 0:.1:10;
% uu2     = 1+.1*randn(length(tt2),1);
% [yg,tg] = lsim(G,uu2,tt2);
% [ys,ts] = lsim(Hr,uu2,tt2);
% 
% W = logspace(-2,2,200);
% figure, 
% subplot(221); hold on, axis tight
% plot(tg,yg,'-')
% plot(ts,ys,'--')
% xlabel('Time [s]'); ylabel('Amplitude [.]');
% legend({'Original model','Identified model'},'Location','Best')
% subplot(223); hold on, axis tight
% bodemag(G,'-',Hr,'r--')
% legend({'Original model','Identified model'},'Location','Best')
% %mdspack.bodemag(G,W,'LineWidth',3)
% %mdspack.bodemag(Hr,W,'--','LineWidth',3)
% subplot(2,2,[2 4]); hold on, grid on
% eigG    = eig(G);
% eigHr   = eig(Hr);
% plot(real(eigG),imag(eigG),'o')
% plot(real(eigHr),imag(eigHr),'.')
% xlabel('Real'); ylabel('Imag.');
% legend({'Original model','Identified model'},'Location','Best')


