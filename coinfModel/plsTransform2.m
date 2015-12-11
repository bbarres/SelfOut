function [Stf] = plsTransform2(S,lam,W,mju,sd)

nS = size(S,2);
nsim = size(S,1);
ins = lam>0;
lamr = repmat(lam(ins),nsim,1);
Scx = S;
Scx(:,ins) = (S(:,ins).^lamr-1)./lamr;
Scx = (Scx - repmat(mju,nsim,1))./repmat(sd,nsim,1);
ncomp = size(W,2);
Stf = zeros(nsim,ncomp);

for i = 1:nsim
    Stf(i,:) = Scx(i,:)*W;
end

