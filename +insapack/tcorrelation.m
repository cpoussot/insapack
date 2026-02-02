% <strong>INSAPACK</strong>
% Main interface for cross correslation
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
%
% Description
% Computes a cross correlation between u and y signals.
%
% Syntax
%  Rxy = insapack.tcorrelation(x,y,KEY)
% 
% Input arguments
%  - x : signal #1 (vector)
%  - y : signal #2 (vector)
%
% Optional arguments
%  - KEY : select from 'biased' of 'un-biased' mode
%
% Output arguments
%  - Rxy : correlated output
% 

function Rxy = tcorrelation(x,y,KEY)

if nargin < 2
    y = x;
end
if nargin < 3
    KEY = 'B';
end
% x = x - mean(x);
% y = y - mean(y);
N = length(x);

switch KEY
    case 'B' % biased 
        alpha = @(m) 1/N;
    case 'U' % un-biased
        alpha = @(m) 1/(N-m);
end
%
Rxy = zeros(N,1);
for m = 1:N
    Rxy(m) = alpha(m-1) * x(1,1:N-m)*y(1,m+1:N)';
    %Rxy(m) = alpha(m-1) * x(1,1:N-m)'*y(1,m:N);
end
