function [f,U0,Y0,G0,sigU2,sigY2,sigUY2,sigG2,b,rho] = non_param_freq(u,y,Ts)

%%% Number of repeated periods of experiments
[nt,Pu] = size(u);
[~,Py]  = size(y);
P       = max(Pu,Py); 
% Treat the case where u or y is a vectors only
if Pu == 1
    u = repmat(u,[1 P]);
end
if Py == 1
    y = repmat(y,[1 P]);
end

%%% FFT
f   = linspace(0,1,nt).'/Ts;
U   = zeros(nt,P);
Y   = U;
for i = 1:P
    U(:,i) = fft(u(:,i))/nt;
    Y(:,i) = fft(y(:,i))/nt;
end

%%% Stochastic analysis of periodic excitations
% >> Average
U0      = mean(U,2);
Y0      = mean(Y,2);
% >> Variance
sigU2   = 1/(P-1)*sum(abs(U-U0).^2,2); % sig^2 = var(y)
sigY2   = 1/(P-1)*sum(abs(Y-Y0).^2,2); % sig^2 = var(u)
% >> Covariance
sigUY2  = 1/(P-1)*sum((Y-Y0).*conj(U-U0),2); 
% >> Correlation: measures the linear relation between two noise sources
rho     = sigUY2./(sqrt(sigU2).*sqrt(sigY2));
% >> FR
% Smoothing the empirical transfer function estimate over successive realizations
%G0  = Y0./U0; % Noob mode
G0  = sum((Y.*conj(U)),2)./sum(U.*conj(U),2);
% >> Bias estimation
b   = -G0.*exp(-P*abs(U0).^2./sigU2) .* ... 
      ( 1-rho.*(U0./sqrt(sigU2)) ./ (Y0./sqrt(sigY2)) );
% >> Variance: the variability of G(k) around its expected value is 
% characterized by the variance sigG2. The larger, the wider the spread 
% around the expected value. For a Gaussian pdf, the interval:
% [-1.96sig;+1.96sig] corresponds to the 95% confidence interval.
sigG2   = 1/P*abs(G0).^2.*( sigY2./abs(Y0).^2 + ...
                            sigU2./abs(U0).^2 - ...
                            2*real(sigUY2./(Y0.*conj(U0))) );
