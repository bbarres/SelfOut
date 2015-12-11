load coinfNR2 
thetaHR12 = thetaHR500;
load coinf13NR1 thetaHR100 thetaHR500
thetaHR13 = thetaHR500;

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
    fR12 = ksdensity(thetaHR12(:,i),x);
    fR13 = ksdensity(thetaHR13(:,i),x);
    
    
    subplot(nrow,ncol,i);
    hold on
    plot(x,fAll2);
    plot(x,fR12,'r');
    plot(x,fR13,'g');

    yla = 1.02*max([fAll2 fR12 fR13]);
    
       
    axis([min(x) max(x) 0 yla]);
    if i==ncol
        legend('Prior','2012','2013')
    end
    
    xlabel(parametrit{i});
    
end