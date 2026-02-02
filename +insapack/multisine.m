% <strong>INSAPACK</strong>
% Main interface for multi-sine signal generation
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
%
% Description
% Computes a multi-sine signal sweeping over FBND frequency range.
%
% Syntax
%  [uc,tc,info] = insapack.multisine(Ns,Ts,FBND,RPHI,ODD,REV,SHOW)
% 
% Input arguments
%  - Ns   : number of samples (integer)
%  - Ts   : sampling period [s] (>0 real)
%  - FBND : frequency range [Hz] (2x1 real vector)
%  - RPHI : random phase (boolean)
%           - false: Schroeder multisine
%           - true: random multisine
%  - ODD  : frequency excited (used for nonlinear detection)
%           - 'all': all frequencies excited
%           - 'odd': odd frequencies excited (even discarded)
%           - 'odd-odd': odd-odd frequencies excited
%           - 'odd-rnd': odd frequencies excited (random block not)
%  - REV  : reverse signal (boolean)
%           - false: fmin to max 
%           - true: fmax to fmin
%  - SHOW : plot signal (boolean)
% 
% Output arguments
%  - uc   : chirp signal in [-1,1]
%  - tc   : time samples
%  - info : additional information
% 

function [uc,tc,info] = multisine(Ns,Ts,FBND,RPHI,ODD,REV,SHOW)

% Nyquist check
Fs  = 1/Ts;
if max(FBND) >= Fs/2
    error('Nyquist frequency issue: check "max(fBnd) < 1/(2Ts)".')
end
%
tc  = (0:Ns-1)*Ts;
f0  = Fs/Ns;
% Select bounds
fId = 1 + round(FBND./f0);
%F   = length(fId(1):fId(2));
F   = fId(2) - fId(1) + 1;
amp = zeros(1,fId(2));
amp(end-F+1:end) = 1;
N       = length(amp);
mask    = zeros(1,length(amp));
%
switch lower(ODD)
    case 'all'
        mask = ~mask;
    case 'odd'
        for k = 1:N/2
            mask(2*k-1) = 1;
        end
        amp = mask.*amp;
    case 'odd-odd'
        for k = 1:N/4
            mask(4*k-3) = 1;
        end
        amp = mask.*amp;
    case 'odd-rnd'
        for k = 1:N/2
            mask(2*k-1) = 1;
            if mod(k,3) == 0
                if rand(1) > .2
                    mask(2*k-1) = 0;
                end
            end
        end
        amp = mask.*amp;
    otherwise
        
end
% 
u   = zeros(F,Ns);

for k = 1:length(amp)
    if RPHI
        phi_k   = rand(1)*2*pi;  % Random multisine
    else
        phi_k   = -k*(k-1)*pi/F; % Schroeder multisine
    end
    if REV
        u(k,:)  = amp(k)*1/F*cos(2*pi*f0*k*tc(end:-1:1) + phi_k);
    else
        u(k,:)  = amp(k)*1/F*cos(2*pi*f0*k*tc           + phi_k);
    end
end
uc      = sum(u,1);
uc      = uc/max(uc(:));
% FFT
L       = numel(uc);
FTuc    = fft(uc)/L;
fc      = linspace(0,1,L)*Fs;
% Info
info.Ts     = Ts;
info.Fs     = Fs;
info.fc     = fc;
info.FTuc   = FTuc;
%
info.uk     = u;
info.mask   = mask;
% Plot
if SHOW 
    FONT_SZ     = 16;
    FONT_SZ2    = 14;
    %
    figure, 
    subplot(211); hold on, grid on, axis tight
    plot(tc,uc,'-','LineWidth',3),
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$t$ [s]','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{u}(t)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Continuous-time'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    %
    subplot(212); hold on; grid on, axis tight
    plot(fc,abs(FTuc),'LineWidth',3),
    hh = gca;
    plot([1 1]*Fs/2,[hh.YLim(1) hh.YLim(2)],'k:','LineWidth',3), 
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$f$ [Hz]','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{U}(f)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Continuous-time (FFT)','Nyquist frequency'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    %
    sgtitle(['Multisine signal $\{N_s,T_s,T_f\}=\{' num2str(Ns) ',' num2str(Ts)  ',' num2str(tc(end)) '\}$'],'Interpreter','latex','Fontsize',20)
end
