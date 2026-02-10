function [la,mu,W,V,R,L] = data2loewner(puls,G0)

id1     = 2:2:length(puls);
id2     = 1:2:length(puls);
%
la_     = 1i*puls(id1);
mu_     = 1i*puls(id2);
Gla_    = G0(id1);
Gmu_    = G0(id2);
%
la = []; W = [];
for ii = 1:length(la_)
    la = [la la_(ii)  conj(la_(ii))];
    W  = [W  Gla_(ii) conj(Gla_(ii))];
end
mu = []; V = [];
for ii = 1:length(mu_)
    mu = [mu mu_(ii)  conj(mu_(ii))];
    V  = [V  Gmu_(ii) conj(Gmu_(ii))];
end
WW(1,1,:) = W; W = WW;
VV(1,1,:) = V; V = VV;
%
k   = length(la);
q   = length(mu);
R   = ones(1,k);
L   = ones(q,1);
%