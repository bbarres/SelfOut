function S = laskeSummaryt4(G,X,RA,AA);

nx = size(X,2);

S = zeros(21+5*nx,1);

ylis = 1e3;

np = size(G,1);
npl = size(G,2);
ninf = [sum(G>-1,2) repmat(npl,np,1)];
if any(ninf(:,1)>1 & ninf(:,1)<npl)
    [B,~,stats] = glmfit(X,ninf,'binomial','constant','off');
    if all(abs(B)<ylis)
        S(1:nx) = B;
        S(nx+1) = std(stats.resid);
    end
end
kohta = nx+1;
multi = [sum(G==0,2) repmat(npl,np,1)];
if any(multi(:,1)>1 & multi(:,1)<npl)
    [B,~,stats] = glmfit(X,multi,'binomial','constant','off');
    if all(abs(B)<ylis)
        S(kohta+(1:nx)) = B;
        S(kohta+nx+1) = std(stats.resid);
    end
end

kohta = kohta+nx+1;
nstrains = zeros(np,1);
for k = 1:np;
    rivi = G(k,:);
    strainsk = unique(rivi(rivi>0));
    nstrains(k) = max(length(strainsk),2*any(rivi==0));
end
if any(nstrains(:,1)>1 & nstrains(:,1)<npl)
    [B,~,stats] = glmfit(X,nstrains,'poisson','constant','off');
    if all(abs(B)<ylis)
        S(kohta+(1:nx)) = B;
        S(kohta+nx+1) = std(stats.resid);
    end
end
kohta = kohta+nx+1;

ninfp = sum(G>-1,2)>0;
S(kohta+1) = mean(ninfp);
S(kohta+2) = std(ninfp);
kohta = kohta+2;

rains = RA>-100;

if sum(rains)>0
    S(kohta+1) = mean(RA(rains));
    if sum(kohta+2)>2
        S(19) = std(RA(rains));
        
        if length(unique(sum(X(rains,:),2)))>4
            [B,~,R] = regress(RA(rains),X(rains,:));
            if all(abs(B)<ylis)
                S(kohta+2+(1:nx)) = B;
                S(kohta+nx+3) = sqrt(mean(R.^2));
            end
        end
    end
end

kohta = kohta+nx+3;

aains = RA>-100;

if sum(aains)>0
    S(kohta +1) = mean(AA(aains));
    if sum(aains)>2
        S(kohta+2) = std(AA(rains));
        
        if length(unique(sum(X(aains,:),2)))>4
            [B,~,R] = regress(AA(aains),X(rains,:));
            if all(abs(B)<ylis)
                S(kohta+2+(1:nx)) = B;
                S(kohta+nx+3) = sqrt(mean(R.^2));
            end
        end
    end
end

kohta = kohta + nx + 3;

for i = 1:4
    S(kohta+i) = sum(RA==i);
    S(kohta +4+i) = sum(AA==i);
end
kohta = kohta + 8;

coinf = (sum(G==0,2)>0);
S(kohta +1) = sum(coinf)/np;
if sum(ninfp)>0
    S(kohta +2) = sum(coinf)/sum(ninfp);
end

S = S';