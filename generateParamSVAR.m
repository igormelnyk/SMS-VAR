%function to generate switching VAR model parameters for simulations
function Dist = generateParamSVAR(Nx, Nm, Ny, Lseq)


numModes = 2^Nm;

Dist.numModes = numModes;

%===================== mode duration distribution =========================
muArray = (Lseq/3)*rand(1,numModes);
Dist.modeDurDist = ModeDurDist(muArray);


%===================== mode transition distribution =======================
sparsity = 0.1;
Tm = sprand(numModes, numModes, sparsity);

%make sure diagonal is zero to prevent self-transitions
ind = 1:(numModes+1):(numModes*numModes);
Tm(ind) = 0;

S = sum(Tm);

%find columns of all-zeros and generate some non-zero colum
ind_zeros = find(S==0);

for i=1:length(ind_zeros)
    
    while(1)
        randCol = sprand(numModes, 1, sparsity);
        randCol(ind_zeros(i))=0;
        if nnz(randCol) > 0
            break;
        end
    end
    
    Tm(:, ind_zeros(i)) = randCol;
end

%S(ind_zeros) = 1;
%ind = sub2ind(size(Tm), ind_zeros, ind_zeros);
%Tm(ind) = 1;

%do summation again (no zeros should appear now)
S = sum(Tm);
assert(isempty(find(S==0))==1);

D = spdiags((1./S)', 0, numModes, numModes);

%normalize
Tm = Tm*D;
Dist.modeTransDist = ModeTransDist(Tm);


%====================== phase transition distribution =====================
Tp = full(sprand(Nx*Nx*numModes, 1, 0.8));
Tp = reshape(Tp, [Nx, Nx, numModes]);

%remove all-zeros distributions
S = squeeze(sum(Tp,1));
ind_zeros = find(S==0);

for i=1:length(ind_zeros)
   
   [f1, f2] = ind2sub(size(S), ind_zeros(i));
   Tp(f1, f1, f2) = 1;
   
end

S = sum(Tp,1);
Tp = Tp./repmat(S, [Nx, 1, ones(1, Nm)]);

%USE T IN LOG SPACE for phase transition distribution
%so that we don't need to transform it to log space multiple times
Dist.phaseTransDist = PhaseTransDist(log(Tp), 1);


%==================== observation transition distribution =================
As = zeros(Ny, Ny, Nx);

%make A stable    
for i=1:Nx    
    
    %A = rand(Ny, Ny);
    A = full(sprandn(Ny, Ny, 0.5));
    A = 0.5*(A+A');
    [U S] = eig(A);
    
    s = diag(S);
    s = s + abs(s(1)) + 0.1;
    s = s./(max(s) + 0.5);
    
    tmp = U*diag(s)*U';
    tmp(abs(tmp) < 1e-5) = 0;
    tmp = 0.5*(tmp' + tmp);
    
    As(:,:,i) = tmp;
end


[V, D] = eig(As(:,:,1));
d = diag(D);
d(1) = 0.99;
ind = find(d(2:end) < 0.2)+1;
d(ind) = 0.6;
D = diag(d);
v = V(:,1);
A = V*D*V'; As(:,:,1) = 0.5*(A'+A);

O = [v rand(Ny,Ny-1)];
[Q, ~] = qr(O);
d = rand(Ny,1);
d(1) = 0;
ind = find(d(2:end) < 0.2)+1;
d(ind) = 0.6;
D = diag(d);
A = Q*D*Q'; As(:,:,2) = 0.5*(A'+A);


O = [v rand(Ny,Ny-1)];
[Q, ~] = qr(O);
d = rand(Ny,1);
d(1) = 0.99; 
ind = find(d(2:end) < 0.2)+1;
d(ind) = 0.6;
D = diag(d);
A = Q*D*Q'; As(:,:,3) = 0.5*(A'+A);

Dist.obsTransDist = ObsTransDist(As);
Dist.v = v;

















