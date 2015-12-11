figure;
plot([0.5 11.5],[1 1],'k:');
hold on

parametrit = {'\gamma','log(\tau^2)','\beta_{1}','\beta_{2}','\beta_{3}',...
    '\alpha_{M,1}','\alpha_{M,2}','\alpha_{M,3}','\alpha_{N,1}', ...
    '\alpha_{N,2}','log(\alpha_{N,3})'};
rrR500Ta = zeros(nitertest,nth);
for i = 1:nth;
    thetaR500i = squeeze(thetaR500(:,i,:));
    thobsi = thetarawobs(:,i)';
    
    eroR500i = thetaR500i-repmat(thobsi,500,1);
    eroPr = repmat(thetaraw(:,i),1,nitertest)-repmat(thobsi,niter,1);
    
    rmseR500i = sqrt(mean(eroR500i.^2));
    rmsePr = sqrt(mean(eroPr.^2));
    
    rrR500Ta(:,i) = rmseR500i./rmsePr;
end

boxplot(rrR500Ta,'labels',parametrit);    
xlabel('Parameter')
ylabel('Relative RMSE');
axis([0.5 11.5 -0 1.25]);

set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1 1.2])
set(gca,'TickLabelInterpreter','tex')