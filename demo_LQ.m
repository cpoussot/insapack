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
% The script uses the "+insapack" matlab package.
% 
clear all, close all, clc;
set(groot,'DefaultFigurePosition', [300 100 1000 600]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for iloe1 = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(iloe1)},'factory','default'),'latex'); end

%%% INSAPACK
addpath('/Users/charles/Documents/GIT/insapack')

%%% >> Chose a SISO model (2 cases are proposed)
nom = 'tf4'
%nom = 'random'
switch lower(nom)
    case 'tf4'
        Gs    = tf([1 9 6],conv([1/100 2*.3/10 1],[10 5 1]));
        Gs_ss = ss(Gs);
    case 'random'
        Gs_ss = stabsep(rss(30));
        Gs = tf(Gs_ss);
end

%%% >> LQ control with infinite horizon, u = -K*(x-x_ref)
% Design parameters
[A,B,C,D]   = ssdata(Gs_ss);
n           = length(A);
nu          = size(B,2);
R           = eye(nu); 
%
rho     = 1e3; % You may play with rho>0
Q       = rho*eye(n);
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
theta = linspace(0,2*pi,1e3);
plot(cos(theta),sin(theta),'g-'); 
plot(cos(theta)-.5,sin(theta),'r--'); 
legend('$L(\imath \omega)$','PM = $L(\imath \omega)$ crosses circle','$L(\imath \omega)$ outside circle','Location','Best'), axis equal
% %% Step response
% figure, hold on
% step(Gs,'b-')
% step(CLs,'r-')