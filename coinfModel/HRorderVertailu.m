

load coinfNR 
thetaHR1 = thetaHR500;
load coinfNRC thetaHR500
thetaHR2 = thetaHR500;
load coinfNRD thetaHR500
thetaHR3 = thetaHR500;
load coinfNRM100 thetaHR500
thetaHR4 = thetaHR500;


parametrit = {'\gamma','log(\tau^2)','\beta_{1}','\beta_{2}','\beta_{3}',...
    '\alpha_{M,1}','\alpha_{M,2}','\alpha_{M,3}','\alpha_{N,1}', ...
    '\alpha_{N,2}','log(\alpha_{N,3})'};
nrow = 3;
ncol = 4;

figure
for i = 1:nth;
    atheta = linspace(thmin(i),thmax(i),1e5);
    atheta([1 end]) = [];
    [fAll,x] = ksdensity(atheta);
    sapu = max(fAll);
    fAll2 = zeros(size(fAll));
    inxit = (x>thmin(i))&(x<thmax(i));
    fAll2(inxit) = sapu;
    fR1 = ksdensity(thetaHR1(:,i),x);
    fR2 = ksdensity(thetaHR2(:,i),x);
    fR3 = ksdensity(thetaHR3(:,i),x);
    fR4 = ksdensity(thetaHR4(:,i),x);
    
    
    subplot(nrow,ncol,i);
    hold on
    plot(x,fAll2);
    plot(x,fR1,'r');
    plot(x,fR2,'g');
    plot(x,fR3,'k');
    plot(x,fR4,'m');
    
    yla = 1.02*max([fAll2 fR1 fR2 fR3]);
    
       
    axis([min(x) max(x) 0 yla]);
    if i==ncol
        legend('Prior','Random','Prevalence','Weighted','M100')
    end
    
    xlabel(parametrit{i});
    
end