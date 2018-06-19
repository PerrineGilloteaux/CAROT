clear all, close all;
PO2 = [31 18 4 2]; %Partial pressure of O2 from blood wall
Pd = [ 5 22 40 69]; % Corresponding distance in um
PO2 = PO2/760*1; % PUt pressure in % O2
x = Pd(1:end)/15; % to fit pixels
y = PO2(1:end);
myfittype = fittype(' scale*exp(-(x.^2)/(2*sigma*sigma))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'scale','sigma'})
myfit = fit(x',y',myfittype,'StartPoint', [ 0.025,2] )
plot(myfit,x,y)

