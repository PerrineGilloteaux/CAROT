function [soldosetotal ] = plotfigure5sigmoid( fracD,delT,survival )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%define sigmoid
myfittype = fittype('1./(1 + exp(-lambda*x + shift))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'lambda','shift'})
index=1;
soldosetotal=nan(4,1);
x= unique([ 1 round(20/fracD) round(linspace(21/fracD,80/fracD,30))])*fracD;
% plot( x,  percent1 ); hold on;
% plot( x,  percent2 ); hold on;
% plot( x,  percent3 ); hold on;
% plot(x,  percent4 ); hold on;
syms dosetotal
for numD = unique([ 1 round(20/fracD) round(linspace(21/fracD,80/fracD,30))])
    load(['Results\cell_num_Dose3D_f200_cl2_m1_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat']);
    load(['Results\cell_num_Dose3D_f200_cl2_m2_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat']);
    load(['Results\cell_num_Dose3D_f200_cl2_m3_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat']);
    load(['Results\cell_num_Dose3D_f200_cl2_m4_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat']);
    percent1(index)=mean(tc1);
    percent2(index)=mean(tc2);
     percent3(index)=mean(tc3);
      percent4(index)=mean(tc4);
    index=index+1;
end
legend_str = {'Perf+HYPO+ECDeath', 'Perf+ Hypo','Hypo alone','no O2 effect'};
for m=1:4
subplot(2,2,m);
percent = matlab.lang.makeValidName(['percent',num2str(m)]);
percent=eval(percent);
[myfit,gof] = fit(x', percent',myfittype,'StartPoint', [1,0] );
if gof.rsquare>0
    eqn = ( 1./(1 + exp(-myfit.lambda*dosetotal + myfit.shift)) == survival);
    soldosetotal(m) = solve(eqn,dosetotal);
    plot(myfit,x,percent)
     title(['Model ',legend_str{m}]);
else
    plot(x,percent,'.b');
    title(['Model ',legend_str{m}, ' fit has failed']);
end
axis([0 90 0 1]);
% 
% title(['dose ',num2str(fracD),'Gy with delT=',num2str(delT)]);
% legend_str = {'Perf+HYPO+ECDeath', 'Perf+ Hypo','Hypo alone','no O2 effect'};
% legend(legend_str);
 xlabel('Total Dose in Gy');
 ylabel('Tumor control in % on 10 samples only');

end

