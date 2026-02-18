% <strong>INSAPACK</strong>
% Main interface for MLBS signal generation
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
%
% Description
% Computes a MLBS signal over FBND frequency range.
%
% Syntax
%  [uc,tc,info] = insapack.mlbs(Ns,Ts,FBND,REV,SHOW)
% 
% Input arguments
%  - Ns   : number of samples - estimated (integer)
%           Ns may be modified to ensure the number of MLBS cell generating 
%           is a 2^D-1 sequence, where D = round(log2(N*Ts/Tc))
%  - Ts   : sampling period [s] (>0 real)
%  - FBND : frequency range [Hz] (2x1 real vector)
%  - REV  : reverse signal (boolean)
%           - false: fmin to max 
%           - true: fmax to fmin
%  - SHOW : plot signal (boolean)
% 
% Output arguments
%  - uc   : chirp signal in {-1,1}
%  - tc   : time samples
%  - info : additional information
% 
% See also for the Polynomial feedback sequence and explanations on PRBS
% and white noise generation:
% https://www.digikey.fr/fr/articles/techzone/2018/mar/use-readily-available-components-generate-binary-sequences-white-noise
% 

function [uc,tc,info] = mlbs(Ns,Ts,FBND,REV,SHOW)

% Nyquist check
Fs  = 1/Ts;
if max(FBND) >= Fs/2
    error('Nyquist frequency issue: check "max(fBnd) < 1/(2Ts)".')
end
% MLBS can only handle frequency range [0 fmax],  
if min(FBND) > 0
    warning('MLBS feature can only deal with frequency bound [0 fmax]. Consider using multi-sine or chirp functions when minimal frequency bound is greater than 0.')
end
%
Tc  = 1/max(FBND);
Nc  = Ns*Ts/Tc;
if Nc < 2
    error('Not enougth points: (i) consider larger "Ns", "Ts" or (ii) smaller "fmax".')
end
% Plot
if nargin < 4
    SHOW = [];
end
if isempty(SHOW)
    SHOW = false;
end
% Generation of simulation time
Fc  = 1/Tc;
D   = round(log2(Nc));
Nc  = 2^D-1;
td  = (0:Nc-1)*Tc;
n   = length(td);
% Generate PRBS length and feedback polynomial function
D   = max(min(D,32),1);
switch D
    case 2
        F = [2 1]; 
    case 3
        F = [3 2]; 
    case 4
        F = [4 3];
    case 5
        F = [5 3]; 
    case 6
        F = [6 5]; 
    case 7
        F = [7 6]; 
    case 8
        F = [8 6 5 4]; 
    case 9
        F = [9 5];
    case 10
        F = [10 7]; 
    case 11
        F = [11 9];
    case 12
        F = [12 6 4 1];
    case 13
        F = [13 4 3 1];
    case 14
        F = [14 5 3 1]; 
    case 15
        F = [15 14]; 
    case 16
        F = [16 15 13 4]; 
    case 17
        F = [17 14]; 
    case 18
        F = [18 11]; 
    case 19
        F = [19 6 2 1];
    case 20
        F = [20, 17]; 
    case 21
        F = [21, 19];
    case 22
        F = [22, 21];
    case 23
        F = [23, 18];
    case 24
        F = [24, 23, 22, 17];
    case 25
        F = [25, 22];
    case 26
        F = [26, 6, 2, 1];
    case 27
        F = [27, 5, 2, 1];
    case 28
        F = [28, 25];
    case 29
        F = [29, 27];
    case 30
        F = [30, 6, 4, 1];
    case 31
        F = [31, 28];
    case 32
        F = [32, 22, 2, 1];
end
% Initial sequence
%si = ones(1,D)
si = randi(2,1,D)-1;
if sum(si) == 0
    si = ones(1,D);
end
% Iteration using feedback polynomial (periodicity 2^(Ns)-1)
ud  = zeros(1,n);
for ii = 1:n
    sN = xor(si(F(1)),si(F(2)));
    for jj = 3:length(F)
        sN = xor(sN,si(F(jj)));
    end
    si       = [sN si(1:D-1)];
    ud(ii)   = si(D);
end
% {-1,+1}
ud = ud(1:n)*2-1;
% Discrete FFT
L       = numel(ud);
FTud    = fft(ud)/L;
fd      = linspace(0,1,L)*Fc;
% Continuous up-sampled & FFT
ups     = Fs;
tc      = 0:Ts:(Nc*Tc-Ts);
L       = numel(tc);
uc      = zeros(1,L);
for ii = 1:L
    if ceil(ii/ups/Tc) > Nc
        uc(ii) = uc(ii-1);
    else
        uc(ii) = ud(ceil(ii/ups/Tc));
    end
end
if REV
    uc = fliplr(uc);
end
FTuc    = fft(uc)/L;
fc      = linspace(0,1,L)*Fs;
% Autocorrelation 
Ruu     = insapack.tcorrelation(uc,uc);
% Output informations
info.td     = td;
info.ud     = ud;
info.fd     = fd;
info.FTud   = FTud;
info.fc     = fc;
info.FTuc   = FTuc;
%
info.Ruu    = Ruu;
info.D      = D;
% 
if SHOW
    FONT_SZ     = 20;
    FONT_SZ2    = 18;
    %
    figure, 
    col = colororder;
    subplot(221); hold on, grid on, axis tight
    plot(td,ud,'o','LineWidth',3),
    plot(tc,uc,'-','LineWidth',3), 
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$t$','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{u}(t)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Discrete-time','Continuous-time (ZOH)'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    % shift the frequency 0 
    zoh     = abs(FTud.*sin(pi*fd/Fc)./(pi*fd/Fc));
    zoh(1)  = FTud(1);
    subplot(222); hold on; grid on, axis tight
    stem(fd,abs(FTud),'LineWidth',3,'Color',col(1,:)), 
    stem(fd,zoh,'LineWidth',3,'Color',col(3,:),'Marker','+')
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$f$','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{U}(f)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Discrete-time FFT','Discrete-time (ZOH) FFT'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    hh = gca;
    %
    subplot(223); hold on; grid on, axis tight
    plot(tc,Ruu,'-','LineWidth',3,'Color',col(2,:)),
    plot(tc,-1/D*ones(1,length(tc)),'k:','LineWidth',3)
    plot([1 1]*Tc,[-1/D 1],'k--','LineWidth',3)
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$\tau$','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{R}_{uu}(\tau)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Continuous-time autocorr.','$-1/N$','$T_c$'},'Location','NorthEast','Interpreter','latex','FontSize',FONT_SZ)
    %
    K       = floor(Tc/Ts);
    xzer    = (Fc:Fc:K*Fc);
    subplot(224); hold on; grid on, axis tight
    plot(fc,abs(FTuc),'-','LineWidth',3,'Color',col(2,:))
    plot(xzer,xzer*0,'.','LineWidth',3,'Color','k','MarkerSize',20)
    set(gca,'XScale','log','YLim',hh.YLim*1.1,'TickLabelInterpreter','latex','FontSize',FONT_SZ)
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$f$','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{U}(f)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Continuous-time FFT','Zeros (mult. $kF_c$)'},'Location','Best','Interpreter','latex','FontSize',FONT_SZ)
    %
    sgtitle(['MLBS signal $\{D,N_c,T_c,T_s,T_f\}=\{' num2str(D) ',' num2str(Nc) ',' num2str(Tc)  ',' num2str(Ts)  ',' num2str(tc(end)) '\}$'],'Interpreter','latex','Fontsize',20)
end
