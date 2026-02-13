% Loewner algorithm (tangential version)
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [hr,info] = lf.loewner_tng(la,mu,W,V,R,L,opt)
%  
% Input arguments
%  - la : interpolation points (k x 1, complex)
%  - mu : interpolation points (q x 1, complex)
%  - W  : k-dimensional ny x nu structure of the function/data 
%         evaluated at points "la" (ny x nu x k, complex)
%         W = H(la), where H(.) is the underlying model
%  - V  : q-dimensional ny x nu structure of the function/data  
%         evaluated at points "mu" (ny x nu x q, complex)
%         V = H(mu), where H(.) is the underlying model
%  - R  : right tangential directions (ny x k, complex)
%  - L  : left tangential directions (q x nu, complex)
%  - opt: optional arguments
%    * target: rational order (if integer >=1), 
%              SVD tolerance (if real <1)
%              automatic (if Inf or [])
%    * D     : D-term (ny x nu, complex)
%              /!\ use with discretion
%    * real  : boolean to try (if possible) real realization 
% 
% Output arguments
%  - hr   : approximation model (handle function)
%  - info : structure with informations about the Loewner world
%    * r    : rational order (integer)
%    * nu   : McMilan degree (integer)
%    * isCC : is complex conjugated (boolean)
%             if true, most of the matrices are converted to real
%    * J    : complex to real transformation matrix (k x q, complex)
%    * la   : lambda's column interpolation points (k x 1, complex)
%    * mu   : mu's row interpolation points (q x 1, complex)
%    * sv   : normalized singular values of [LL SS] (min(q,k), real)
%    * sv_nu: normalized singular values of LL (min(q,k), real)
%    * LL   : Loewner matrix (q x k, complex)
%    * SS   : shifted Loewner matrix (q x k, complex)
%    * LA   : Lambda matrix (k x k, complex)
%    * MU   : Mu matrix (q x q, complex)
%    * L,R  : left, right tangential data (data x tangent directions)
%    * V,W  : tangential input and output matrices (q x nu & ny x k, complex)
%    * D    : D-term (ny x nu, complex)
%    * H    : full Loewner form (handle function)
%               H(s)=W(-s LL+SS)\V+D
%    * X,Y  : left and right projectors (k x r & q x r, complex)
%    * Rr,Lr: right and left tangential directions  (nu x r & r x ny, complex)
%    * Hr   : compressed Loewner form (state-space, complex)
%               Hr(s)=Cr(sEr-Ar)\Br+Dr, 
%             where (Er,Ar,Br,Cr,Dr) are available in info.Er ...
%    * lat  : compressed column (right) interpolation points (r x 1, complex)
%    * mut  : compressed row (left) interpolation points (r x 1, complex)
%    * ... Others not documented yet
% 
% Note 
% Sylvester equations 
%    MU*LL-LL*LA = V*R-L*W 
%    MU*SS-SS*LA = MU*V*R-L*W*LA
% may be checked for complex form only, i.e.
% if ~info.isCC
%   test1 = info.MU*info.LL - info.LL*info.LA;
%   test2 = info.V*info.R - info.L*info.W;
%   norm(test1-test2) % small
%   test1 = info.MU*info.SS - info.SS*info.LA;
%   test2 = info.MU*info.V*info.R - info.L*info.W*info.LA;
%   norm(test1-test2) % small
% end
% 
% Description
% Loewner rules.
%

function [hr,info] = loewner_tng(la_,mu_,W_,V_,R,L,opt)

TOL_CC  = 1e-12;
TOL_SV  = 1e-15;
%
if ~isempty(intersect(la_,mu_)) 
    error('Repetition in "la" and "mu"')
end
%
[ny,nu,~]   = size(W_);
if nargin < 7 || ~isa(opt,'struct')
    D           = zeros(ny,nu);
    robj        = inf;
    MAKE_REAL   = true;
elseif isa(opt,'struct')
    if isfield(opt,'target')
        robj = opt.target;
    else
        robj = inf;
    end
    if isfield(opt,'D')
        D = opt.D;
    else
        D = zeros(ny,nu);
    end
    if isfield(opt,'real')
        MAKE_REAL = opt.real;
    else
        MAKE_REAL = true;
    end
end
%
k   = length(la_);
q   = length(mu_);

%%% Reshape data
W   = zeros(ny,k);
V   = zeros(q,nu);
for ii = 1:k
    W(1:ny,ii) = W_(:,:,ii)*R(:,ii);
end
for ii = 1:q
    V(ii,1:nu) = L(ii,:)*V_(:,:,ii);
end

%%% Complex conjugation
isCC = false;
if (abs(sum(imag(la_)))/max(abs(la_))  <TOL_CC) && ...
   (abs(sum(imag(mu_)))/max(abs(mu_))  <TOL_CC) && ...
   (abs(sum(imag(W(:))))/max(abs(W(:)))<TOL_CC) && ...
   (abs(sum(imag(V(:))))/max(abs(V(:)))<TOL_CC) && ...
   (q==k) && ...
   MAKE_REAL
    isCC = true;
end
if MAKE_REAL & ~isCC
    warning on
    warning('/!\ Real form cannot be computed, check la, mu, L and R')
end

%%% Loewner matrices
LL  = zeros(q,k);
SS  = LL;
for ii = 1:q
    for jj = 1:k
        num1        = V(ii,:)*R(:,jj) - L(ii,:)*W(:,jj);
        num2        = mu_(ii)*V(ii,:)*R(:,jj) - L(ii,:)*W(:,jj)*la_(jj);
        den         = mu_(ii)-la_(jj);
        LL(ii,jj)   = num1/den;
        SS(ii,jj)   = num2/den;
    end
end
%norm(-SS-SS') 

%%% D-term
if ~norm(D) == 0
    SS  = (SS - L*D*R);
    V   = V - L*D;
    W   = W - D*R;% - ou + ???
end

%%% Go real
if isCC
    J0  = (1/sqrt(2))*[1 1i; 1 -1i];
    J   = [];
    kk  = 1;
    while length(J) < length(la_)
        if imag(la_(kk)) == 0
            J   = blkdiag(J,1);
            kk  = kk + 1;
        else
            J   = blkdiag(J,J0);
            kk  = kk + 2;
        end
    end
    LL  = real(J'*LL*J);
    SS  = real(J'*SS*J);
    V   = real(J'*V);
    W   = real(W*J);
end

%%% Compressed model
% orders
[L1,S1,~]   = svd([LL,SS],'econ','vector');
[~,S2,R2]   = svd([SS',LL']','econ','vector');
S           = S1;
if numel(S2) < numel(S1)
    S = S2;
end
S_nu        = svd(LL,'econ','vector');
sv          = S/S(1,1);
sv_nu       = S_nu/S_nu(1,1);
if isempty(robj) | isinf(robj)
    r   = sum(sv>TOL_SV);
    nu_ = sum(sv_nu>TOL_SV);
elseif robj < 1
    r   = sum(sv>robj);
    nu_ = sum(sv_nu>robj);
elseif robj >= 1
    r   = robj;
    nu_ = robj;
end
X = []; Y = [];
if robj > 0
    Y   = L1(:,1:r);
    X   = R2(:,1:r);
    % compressed
    Er  = -Y'*LL*X;
    Ar  = -Y'*SS*X;
    Br  = Y'*V;
    Cr  = W*X;
    Dr  = D;
    Lr  = Y'*L;
    Rr  = R*X;
else
    Er  = -LL;
    Ar  = -SS;
    Br  = V;
    Cr  = W;
    Dr  = D;
    Lr  = L;
    Rr  = R;
end

%%% Output
hr  = @(s) Cr*((Er*s-Ar)\Br) + Dr;
h   = @(s) W*((-LL*s+SS)\V) + D;
Hr  = dss(Ar,Br,Cr,Dr,Er);

%%% Compressed IP
[Tla,LAt,~] = eig(Ar+Br*Rr,Er); % right IP
[~,MUt,Tmu] = eig(Ar+Lr*Cr,Er); % left IP
%[Tla,LAt] = eig(inv(Er)*(Ar+Br*Rr));
%[Tmu,MUt] = eig((Ar+Lr*Cr)*inv(Er));
LAt         = diag(LAt);
MUt         = diag(MUt);
LLt         = Tmu*Er*Tla;
SSt         = Tmu*Ar*Tla;
Vt          = Tmu*Br; % "-"?
Wt          = Cr*Tla; % "-"?
Lt          = Tmu*Lr;
Rt          = Rr*Tla;

%%% Information
% Loewner
info.r      = r;
info.nu     = nu_;
info.isCC   = isCC;
if isCC
    info.J  = J;
end
info.la     = la_(:);
info.mu     = mu_(:);
info.sv     = sv;
info.sv_nu  = sv_nu;
info.LL     = LL; 
info.SS     = SS;
info.LA     = diag(info.la);
info.MU     = diag(info.mu);
info.L      = L;
info.R      = R;
info.V      = V;
info.W      = W;
info.D      = D;
info.H      = h;
% Compression
info.X      = X;
info.Y      = Y;
info.Rr     = Rr;
info.Lr     = Lr;
info.Hr     = Hr;
info.Er     = Er;
info.Ar     = Ar;
info.Br     = Br;
info.Cr     = Cr;
info.Dr     = Dr;
% Barycentric/modal form
info.lat    = LAt;
info.mut    = MUt;
if (ny == 1) && (nu == 1) 
    info.vt = Vt./Lt; %-(Tmu*Y'*V)./(Tmu*Y'*L);
    info.wt = Wt./Rt; %-(W*X*Tla)./(R*X*Tla);
end
info.LLt    = LLt; 
info.SSt    = SSt;
info.Vt     = Vt;
info.Wt     = Wt;
