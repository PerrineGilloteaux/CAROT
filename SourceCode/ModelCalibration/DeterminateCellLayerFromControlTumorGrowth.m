clear all, close all;
%% Explain stereo here:
% Tumor disk surface=Pi*R^2~=volcell*ncells; we will approximate as equal
%<=> R=sqrt((volcell*ncell)/pi)
% V=4/3*pi*R^3
% <=> V= 4/3*volcell*sqrt(volcell/pi)*(ncells*sqrt(ncells)
%<=> V= A*ncells *sqrt(ncells)
% we will call A volcell in the following for simplification, as a scale
% factor
%% Start to generate appropriate data
%% Data for growth control only
fracD=0;
numD=0;
delT=1;
num_sample=100;
field_size=400;
RTres=0.01; % At beginning, difonly 1% of tumor cells are hypoxic.
O2_th=0.002; % 0.2% Oxygen level defines hypoxia below it.
for n_cell_layers=1:5


fprintf('%d Gy, %d fractions, %d days between exposure\n',fracD,numD,delT);
cell_num1_arr=[];
leaking_factor=1.5;
waitfor=21;
extradays0=26;
cell_num1_arr=[];
sim_count=0;
parfor i_sample = 1:num_sample; % number of repeat for average
% for i_sample = 1:1;
%       parallelg=0;   
    parallelg=1;
    i_sample
   
    % Tstack not used here for parfor: use a classical for loop to
    % get T_stack1
    [T_stack1,o2_stack1,cell_num1,o2_lev1] = Hypoxia_model_v10(field_size,fracD,numD,delT,RTres,O2_th,leaking_factor,extradays0,0,0,n_cell_layers,parallelg,waitfor);
    
    tc1(i_sample) = sum(cell_num1(end,1:2))==0;
    cell_num1_arr(:,:,i_sample) = cell_num1;
    
    
end
cell_num1 = mean(cell_num1_arr,3);
cell_num1_std = std(cell_num1_arr,[],3);


sim_count = sim_count + 1;
save(['Results\cell_num_Control_ncl_',num2str(n_cell_layers),'D_',num2str(fracD),'_nbdose',num2str(numD),'_delT',num2str(delT),'.mat'],'cell_num1','cell_num1_std','tc1');

end
%% Define Experimental Data reference from Potiron (Figure 1B PlosOne 2013)not fINISHED
%Growth tumor experimental data
data1x=[4
8
9
12
14
16
18
22
];

data1y=[50
75
80
100
140
175
190
230
];
data2x=[0
8
10
15
21
24
28
];

data2y=[350
450
650
900
1250
1500
1900];

data3x=[25
29
33
39
46
54
];
data3y=[80
130
175
250
700
1150
];
data4x=[ 0 14];
data4y=[240 750];
  dt=1;  
plot(dt+data1x-4,data1y,'o');hold on,
plot(dt+data2x-8+30+1,data2y,'o');hold on,
plot(dt+data3x-27+9+1,data3y,'o');hold on,
plot(dt+data4x+21,data4y,'o');hold on



MeasuredDaysirr=[0 1 3 7 14]+dt+21;

TumorWeightirr=[240 220 290 280 240];
plot(MeasuredDaysirr,TumorWeightirr,'-o');

xlabel('Days');
ylabel('Tumor Volume in mm^3');
t=3;
for n_cell_layers=1:5
fracD=0;numD=0;delT=1;
load(['Results\cell_num_Control_ncl_',num2str(n_cell_layers),'D_',num2str(fracD),'_nbdose',num2str(numD),'_delT',num2str(delT),'.mat']);


volcell=48/(sum(cell_num1(t,[1,2,5]),2).*sqrt(sum(cell_num1(t,[1,2,5]),2)));
volcell^(1/3)
errorbar(volcell*sum(cell_num1(:,[1,2,5]),2).*sqrt(sum(cell_num1(:,[1,2,5]),2)), sum(cell_num1_std(:,[1,2,5]),2).*sqrt(sum(cell_num1(:,[1,2,5]),2))*volcell,'r');hold on;
end

