clear all, close all;
PO2 = [21 8 3 1 2.8 3]; %Partial pressure of O2 from blood wall
Pd = [5 50 100 150 200 320]; % Corresponding distance in um
PO2 = PO2/760*1; % PUt pressure in % O2
x = Pd(1:end)/15; % to fit pixels
y = PO2(1:end)/1;
myfittype = fittype(' scale*exp(-(x.^2)/(2*sigma*sigma))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'scale','sigma'})
myfit = fit(x',y',myfittype,'StartPoint', [1 2.4] )
plot(myfit,x,y)

