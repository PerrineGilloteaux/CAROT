clear all;
close all;


soldose=nan(4,10);
proba=0.5;
tabledelT=[  1 1 2 2 3 3 1];
indexD=0;
for fracD= [  2 3 4 6 8 10 1]
    indexD=indexD+1;
    delT=tabledelT(indexD);
    soldose(:,fracD)=plotfigure5sigmoid( fracD,delT,proba );
saveas(gcf,['Figures/tumorcontrol50pcv6_fracD_',num2str(fracD),'_delT',num2str(delT),'.fig']);
saveas(gcf,['Figures/tumorcontrol50pcv6_fracD_',num2str(fracD),'_delT',num2str(delT),'.png']);
end
figure,
legend_str = {'Perf+HYPO+ECDeath', 'Perf+ Hypo','Hypo alone','no O2 effect'};
for m=1:4
    dose=1:10;
   plot(dose(~isnan(soldose(m,:))),soldose(m,~isnan(soldose(m,:))),'o-') ;hold on,
end
legend(legend_str);