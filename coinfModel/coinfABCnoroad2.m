load data2012
XX = Xb(:,1:3);
np = size(XX,1);
XX(:,2:3) = XX(:,2:3) - repmat(mean(XX(:,2:3)),np,1);
niter = 1e5;
nitertest = 1e2;
npatch = size(XX,1);

Sindeksit = [13 35];
Sala = [0.02 0.005];
Syla = [0.7 0.5];
Sabs = 1e2;

thr = 4;

thmin = [-10 -1 -11 -5 -5  -4 0 -2  1 -5 -5 ];
thmax = [10 5 -1 1 3  1 7 11  9 5 5 ];

einx = [2 11];

nth = length(thmin);
nS = 36;

SobsH = laskeSummaryt4(Gobs,XX,Ab(:,2),Ab(:,1));


thetaraw = repmat(thmin,niter,1) + repmat((thmax-thmin),niter,1).*rand(niter,nth);

theta = [thetaraw(:,1:2) ones(niter,1) thetaraw(:,3:end)];
theta(:,einx) = exp(theta(:,einx));
tic
S = zeros(niter,nS);
parfor i = 1:niter
    th = theta(i,:)';
    [GG,RARA,AAAA] = simuloiCoinfection3(th,XX);
    S(i,:) = laskeSummaryt4(GG,XX,RARA,AAAA);
end
toc

ncomp = 20;
testready = false;
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
testready=1;
postProcess4
TestKuva
save coinfNR2


