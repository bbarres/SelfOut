function [Stf,lam,W,mju,sd,pct] = plsTransformB(S,theta,ncomp);

nS = size(S,2);
nsim = size(S,1);

Scx = zeros(size(S));
%lam = zeros(1,nS);
lam = 0.5*ones(1,nS);
mju = zeros(1,nS);
sd = ones(1,nS);
for i = 1:nS;
    %[Scx(:,i),lam(i)] = boxcox(S(:,i));
    if ~any(S(:,i)<0)  
        Scx(:,i) = boxcox(lam(i),S(:,i));
    else
        Scx(:,i) = S(:,i);
        lam(i) = -1;
    end
    mju(i) = mean(Scx(:,i));
    if length(unique(Scx(:,i)))>1
        sd(i) = std(Scx(:,i));
    end
    Scx(:,i) = (Scx(:,i)-mju(i))./sd(i);
    %disp(i);
end

%[~,~,~,~,~,pctvar,mse] = plsregress(Scx,theta,ncomp);
%[~,ncopt] = find(diff(diff(mse(2,:)))>0,1,'First');
ncopt = ncomp;
[~,~,Stf,~,~,pct,~,stats] = plsregress(Scx,zscore(theta),ncopt);
W = stats.W;
%
%  Stf = zeros(nsim,ncopt);
%
%  for i = 1:nsim
%      Stf(i,:) = Scx(i,:)*W;
%  end
