% DEMO LQ CONTROL
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
% Date  : December 2014 (creation)
%         March 2022 (update)
% 
% Description
% Demonstration script for LQ control design and robustness properties.
% This small script illustrates how to construct and analyse the LQ control
% design in a SISO setting. In particular, we compute the LQ controller by
% solving the CARE equation, then we illustrate some classic margin
% properties of this special and well known control method.
% 
% Note
% The script uses the "+insa" matlab package.
% 

clear all, close all, clc

%%% >> Graphic settings
insa.initGraphics

%%% >> Chose a SISO model (2 cases are proposed)
nom = 'tf4'
%nom = 'random'
switch lower(nom)
    case 'tf4'
        Gs    = tf([1 9 6],conv([1/100 2*.1/10 1],[10 5 98 1]));
        Gs_ss = ss(Gs);
    case 'random'
        Gs_ss = stabsep(rss(30));
        Gs = tf(Gs_ss);
end

%%% >> LQ control with infinite horizon, u = -K*(x-x_ref)
% Design parameters
A   = Gs_ss.a;
B   = Gs_ss.b;
C   = Gs_ss.c;
nu  = size(B,2);
n   = length(A);
R   = eye(nu); 
rho = 1e3; % You may play with rho>0
Q   = rho*eye(n);
% Continuous Algebraic Riccati Equation resolution ...
[P,L,G] = care(A,B,Q,R);   % CARE
% ... and controller construction
K       = R\eye(nu)*B'*P;  % Optimal control u = -Kx
% Transfer of interest
Ls  = tf(ss(A,B,K,0));       % Loop transfer => for robustness
Ss  = (1+Ls)\1;              % Sensitivity function => for disturbance rejection
Ts  = 1-Ss;                  % Complementary sensitivity function
h   = 1/(C*((-A+B*K)\B));    % Steady-state closed-loop gain
CLs = h*tf(ss(A-B*K,B,C,0)); % Closed-loop => input-output performance
% Margin and Hinf-norm
gamma = norm(CLs,inf) % Hinf norm
S     = allmargin(Ls) % Margins 
% Frequency responses
figure 
subplot(121), hold on
bodemag(Ss,'b',Ts,'r-',tf(1,1),'k-'), grid on
title('Sensitivity function')
legend('$S(\imath \omega)$','$T(\imath \omega)$','1','Location','Best')
subplot(122), hold on
sigma(CLs,'b',gamma*tf(1,1),'r'), grid on
title('Closed-loop')
legend('$T(\imath \omega)$','$\gamma$','Location','Best')
% Robustness margin illustration
figure
subplot(121), hold on
margin(Ls)
legend('$L(\imath \omega)$','Location','Best')
subplot(122), hold on
nyquist(Ls,'b-'), grid on
insa.circle([ 0,0],1,100,'g-'); 
insa.circle([-1,0],1,100,'r--'); 
legend('$L(\imath \omega)$','PM = $L(\imath \omega)$ crosses circle','$L(\imath \omega)$ outside circle','Location','Best'), axis equal
% %% Step respo
% figure, hold on
% step(Gs,'b-')
% step(CLs,'r-')