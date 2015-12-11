function [G,RA,AA] = simuloiCoinfection3d(param,X)

% theta - sisältää kaikki parametrit.

npatch = size(X,1);
nx = size(X,2);
Mmax = 30;


khi = param(1);
sigZ = param(2);
taumju = param(3);
bet = param(3 + (1:nx));
alfaM = param(3 + nx + (1:nx));
alfaN = param(3 + 2*nx + (1:nx));


sZ = sqrt(sigZ);

%M
em = exp(X*alfaM);
M = poissrnd(em);
M = min(M,Mmax);

en = exp(X*alfaN);
N = poissrnd(en);
Nmax = 1e3;
Nsim = min(N,Nmax);


maxM = max(M);
mju = taumju*randn(npatch,maxM);
npl = 10;
G = -ones(npatch,npl);

nsim = 1e2;
RA = NaN(npatch,1);
AA = NaN(npatch,1);
kerr = 5;

for j = 1:npatch
    Mj = M(j);
    if Mj>0
        
        nj = Nsim(j);
        %L = chol(S(1:Mj,1:Mj),'lower');
        eb =  X(j,:)*bet;
        theta = repmat(eb+mju(j,1:Mj),nj,1);
        Z = theta + sZ*randn(nj,Mj);
        %[~,jarj] = sort(eb+mju(j,1:Mj),'descend');
        
        apug = exp(eb+mju(j,1:Mj));
        apug = apug./sum(apug);
        for l = 1:nj
            gg = gamrnd(kerr*apug,1);
            [~,jarj] = sort(gg,'descend');
            infec = (Z(l,:)>0);
            ekainf = find(infec(jarj),1,'first');
            loput = jarj((ekainf+1):end);
            Z(l,loput) = Z(l,loput) + khi;
        end
        posit = sum(Z>0,2)';
        
        infected = find(posit>0);
        
        ninf = length(infected);
        
        if ninf>=npl
            rivit = infected(1:npl);
        elseif nj<N(j)
            ero = npl-ninf;
            kesken = 1;
            theta2 = repmat(eb+mju(j,1:Mj),ero,1);
            Z2 = theta2;
            pj = 1;
            while kesken
                lisa = sZ*randn(1,Mj);
                if any(lisa+theta2(pj,:)>0)
                    Z2(pj,:) = Z2(pj,:) + lisa;
                    infec = Z2(pj,:)>0;
                    gg = gamrnd(kerr*apug,1);
                    [~,jarj] = sort(gg,'descend');
                    ekainf = find(infec(jarj),1,'first');
                    loput = jarj((ekainf+1):end);
                    Z2(pj,loput) = Z2(pj,loput) + khi;
                    pj = pj+1;
                end
                nj = nj+1;
                if pj>ero || nj >=N(j) || nj>4*Nmax
                    kesken = 0;
                    pj = pj-1;
                end
            end
            Z = [Z; Z2(1:pj,:)];
            posit = sum(Z>0,2)';
            
            rivit = [infected (size(Z,1)+1 - (1:pj))];
            
        else
            rivit = infected;
        end
            
        nplj = length(rivit);
        apu = -ones(1,nplj);
        apu(posit(rivit)>1) = 0;
        ins = find(posit(rivit)==1);
        for k = 1:length(ins);
            in = ins(k);
            apu(in) = find(Z(rivit(in),:)>0);
            
        end
        G(j,1:nplj) = apu;
        osuus = ninf/nj;
        RA(j) = 1+(osuus>0.01) + (osuus>0.25) + (osuus>0.5) + (osuus>0.75);
        maara = osuus*N(j);
        AA(j) = 1+(maara>10) + (maara>100) + (maara>1000);
    end
end

        


