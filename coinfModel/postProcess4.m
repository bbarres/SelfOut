tic


[Stf,lam,W,mju,sd] = plsTransformB(S,thetaraw,ncomp);
if testready
    [Stfobs] = plsTransform2(Sobs,lam,W,mju,sd);
end
[StfobsH] = plsTransform2(SobsH,lam,W,mju,sd);


if testready
    order = zeros(niter,nitertest);
    
    thetaR500 = zeros(500,nth,nitertest);
    
    for i = 1:nitertest;
        Stfobsi = Stfobs(i,:);
        ero = sum((Stf-repmat(Stfobsi,niter,1)).^2,2);
        apu = [ero (1:niter)'];
        apu = sortrows(apu);
        
        thetaR500(:,:,i) = regressPost(thetaraw,Stf,Stfobsi,apu(1:500,2));
        order(:,i) = apu(:,2);
    end
end
ero = sum((Stf-repmat(StfobsH,niter,1)).^2,2);
apu = [ero (1:niter)'];
apu = sortrows(apu);

thetaHR500 = regressPost(thetaraw,Stf,StfobsH,apu(1:500,2));


disp('Post-processing done');
toc