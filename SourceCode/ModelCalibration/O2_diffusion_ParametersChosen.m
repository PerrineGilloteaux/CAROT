% Model of O2 diffusion
% Setting parameters so that O2 profile matches reality
% Pixel size is about 15 um (one cell)
%
clear all, close all;
% gaussf function bugs if multithread
%% Step 1: show O2 diffusion
dipsetpref('NumberOfThreads',1)

% Published profile
% PO2 in air = 160 mmHg, or 21%
% PO2 in vessel = 100 mmHg
% PO2 in tissue: blood wall - 83: 933?963, 2003; 10.1152/physrev.00034.2002.
PO2 = [55 31 18 4 2]; %Partial pressure of O2 from blood wall
Pd = [0 5 22 40 69]; % Corresponding distance in um
PO2 = PO2/760*1; % PUt pressure in % O2

% create image with one vessel in middle of image and apply diffusion


dif_cst = 2.4;
%dif_cst=sqrt(1.806);
Pscale = 1.18;
blood_density=0.037963;
img = newim(15,15);
img(7,7) = 1 ;% Set center at O2 level corresponding to center of vessell

img = min(max(gaussf(img*Pscale,dif_cst),0.001),0.05);
pro = double(img(7,7:end));
plot(0:15:15*7,pro,'r-');hold on ;plot(Pd,PO2,'o'); hold on;
legend('Model used','Experimental Data');
xlabel('Distance in um from vessel centerline (15=vessel wall)');
ylabel('O2 level');
grid on;
title (['Diffusion and Pscale set: demo on one vessel ',num2str(dif_cst)]);



%% Step 2- Check frac_hypo only 1% of tissu is hypoxic (below 0.2% O2) and overal o2 is 1.3%

%dif_cst = 2.1; % diffusion constant used for oxygen. Obtained from another simulation (O2_diffusion.m)
%Pscale = 0.83; % Value at 15 um away to scale on (obtained so that for dif_cst is 2.5, first pixel out is at ~0.027)
field_size = 200;
for i=1:10
rand_array = rand(field_size, field_size);

img = newim(field_size,field_size);
cell_array = zeros(field_size,field_size);
cell_array(rand_array<blood_density)= 1;
cell_img = dip_image(cell_array);
img=min(max(gaussf(cell_img*Pscale,dif_cst),0.001),0.05);
frac_hypo(i) = sum(img<0.002)/field_size^2;
o2_level(i)= mean(img(cell_img==0));
end
title ({['Diffusion and Pscale set: demo on one vessel ',num2str(dif_cst)],['frac hypo: ',num2str(mean(frac_hypo),2),'+/-',num2str(std(frac_hypo),2)],...
    [' mean o2level: ',num2str(mean(o2_level),2),'+/-',num2str(std(o2_level),2)]});


