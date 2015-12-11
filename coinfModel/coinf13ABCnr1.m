%Use same simulations as with 2012 data
load coinfNR2
load data2013;
XX = Xc(:,1:3);
np = size(XX,1);
XX(:,2:3) = XX(:,2:3) - repmat(mean(XX(:,2:3)),np,1);
%niter = 1e5;
nitertest = 1e2;
npatch = size(XX,1);

SobsH = laskeSummaryt4(Gobs,XX,Ac(:,2),Ac(:,1));

ncomp = 20;
postProcess4;

tic
obsw = 2;
thmedian = median(thetaHR500);
thetarawobs = zeros(nitertest,nth);
Sobs = zeros(nitertest,nS);
parfor i = 1:nitertest;
    huono = 1;
    laskuri = 0;
    while huono
        huono = 0;
        thrw =  thmedian - obsw/2 + obsw*rand(1,nth);
        thetarawobs(i,:) = thrw;
        th = [thrw(1:2) 1 thrw(3:end)]';
        th(einx) = exp(th(einx));
        %thetaobs(i,:) = th';
        [GGo,RAo,AAo] = simuloiCoinfection3(th,XX);
        SSS = laskeSummaryt4(GGo,XX,RAo,AAo);
        sin = SSS(Sindeksit);
        if any(sin<Sala | sin>Syla)
            huono =1;
        elseif any(abs(SSS(1:26))>Sabs)
            huono = 1;
        end
        laskuri = laskuri + 1;
        
    end
    disp(laskuri);
    Sobs(i,:) = SSS;
end
toc
disp('Test data done');

testready =1;
postProcess4
save coinf13NR1



