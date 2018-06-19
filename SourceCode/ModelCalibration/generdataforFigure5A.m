% This script will establish tumor control curves as a function total dose
% for various fraction protocol Missing values were 1 gy, 3 gy , 6 gy and
% 10gy
% August 2015
clear all, close all
tic
% %% test 1: Perf+HYPO+ECDeath
% leak_factor = 1.3; %each time there is a death hit in blood vessel, that vessel becomes leak_factor times more leaky
% % Basically, leak_factor the scaling factor for O2
% % point source level at vesssel position
% O2_option=1;
% EC_death_option=1;
% %% test 2: Perf+ Hypo
% leak_factor = 1.3; %each time there is a death hit in blood vessel, that vessel becomes leak_factor times more leaky
% % Basically, leak_factor the scaling factor for O2
% % point source level at vesssel position
% O2_option=1;
% EC_death_option=0;
% %% test 3: Hypo alone
% leak_factor = 1; %each time there is a death hit in blood vessel, that vessel becomes leak_factor times more leaky
% % Basically, leak_factor the scaling factor for O2
% % point source level at vesssel position
% O2_option=1;
% EC_death_option=0;
% %% test 4 : no O2 effect
% leak_factor = 1; %each time there is a death hit in blood vessel, that vessel becomes leak_factor times more leaky
% % Basically, leak_factor the scaling factor for O2
% % point source level at vesssel position
% O2_option=0;
% EC_death_option=0;



field_size = 400; % in v8 tumor dilation is fieldsize/4 not by 2 to avoid space saturation


extra_days = 15; % number of days to track after IR (important to see if cell goes to 0)
RTres=0.01; % At beginning, difonly 1% of tumor cells are hypoxic.
O2_th=0.002; % 0.2% Oxygen level defines hypoxia below it.
num_sample = 10; % number of simulated tumors to be ran to compute average response
%sim_results = zeros(11*4,6); % For each dose and fraction, report fracD, numD, TotalD, Tumor control, TC error, delT
% 4 doses, 11 fraction regimens
tc1 = nan(num_sample,1);
tc2 = nan(num_sample,1);
tc3 = nan(num_sample,1);
tc4 = nan(num_sample,1);

sim_count = 1;
T_stack1=[];
o2_stack1=[];


o2_lev1=[];
indexD=0;
parallel=1;
n_cell_layer=3;
%% Data for growth control only
fracD=0;
numD=0;
delT=1;
fprintf('%d Gy, %d fractions, %d days between exposure\n',fracD,numD,delT);
cell_num1_arr=[];
leaking_factor=1.5;
waitfor=21;
extradays0=26;
parfor i_sample = 1:num_sample; % number of repeat for average
% for i_sample = 1:1;
%       parallelg=0;   
parallelg=1;
    i_sample
   
    % Tstack not used here for parfor: use a classical for loop to
    % get T_stack1
    [T_stack1,o2_stack1,cell_num1,o2_lev1] = Hypoxia_model_v10(field_size,fracD,numD,delT,RTres,O2_th,leaking_factor,extradays0,0,0,n_cell_layer,parallelg,waitfor);
    
    tc1(i_sample) = sum(cell_num1(end,1:2))==0;
    cell_num1_arr(:,:,i_sample) = cell_num1;
    
    % forbidden in parfor save(['Results\cell_num_DoseotherControl3D',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'_run',num2str(i_sample),'.mat']);
end
cell_num1 = mean(cell_num1_arr,3);
cell_num1_std = std(cell_num1_arr,[],3);


sim_count = sim_count + 1;
save(['Results\cell_num_Dose3D_f1920_cl3_m1_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat'],'cell_num1','cell_num1_std','tc1');

    
%% actual simulation for numcell only (generate stack where needed with parralel option set to 0 and a classical for loop
% (one patient is enough for stack by the way)
tabledelT=[  1 ];
for fracD = [  2 ]
    %fracD = [0 ] % Dose per fraction
    indexD=indexD+1;
    delT=tabledelT(indexD);
    for numD = 10 % number of fractions (DOSE TOTAL=70?)
        
        fprintf('%d Gy, %d fractions, %d days between exposure\n',fracD,numD,delT);
        cell_num1_arr=[];
        cell_num2_arr=[];
        cell_num3_arr=[];
        cell_num4_arr=[];
        parfor i_sample = 1:num_sample; % number of repeat for average
            i_sample
            % Tstack not used here for parfor: use a classical for loop to
            % get T_stack1
            [T_stack1,o2_stack1,cell_num1,o2_lev1] = Hypoxia_model_v10(field_size,fracD,numD,delT,RTres,O2_th,leaking_factor,extra_days,1,1,n_cell_layer,parallel,waitfor);
            
            tc1(i_sample) = sum(cell_num1(end,1:2))==0;
            cell_num1_arr(:,:,i_sample) = cell_num1;
            [T_stack1,o2_stack1,cell_num2,o2_lev2] = Hypoxia_model_v10(field_size,fracD,numD,delT,RTres,O2_th,leaking_factor,extra_days,1,0,n_cell_layer,parallel,waitfor);
            
            tc2(i_sample) = sum(cell_num2(end,1:2))==0;
            cell_num2_arr(:,:,i_sample) = cell_num2;
            [T_stack1,o2_stack1,cell_num3,o2_lev1] = Hypoxia_model_v10(field_size,fracD,numD,delT,RTres,O2_th,1,extra_days,1,0,n_cell_layer,parallel,waitfor);
            
            tc3(i_sample) = sum(cell_num3(end,1:2))==0;
            cell_num3_arr(:,:,i_sample) = cell_num3;
            [T_stack1,o2_stack1,cell_num4,o2_lev4] = Hypoxia_model_v10(field_size,fracD,numD,delT,RTres,O2_th,1,extra_days,0,0,n_cell_layer,parallel,waitfor);
            
            tc4(i_sample) = sum(cell_num4(end,1:2))==0;
            cell_num4_arr(:,:,i_sample) = cell_num4;
            % forbidden in parfor save(['Results\cell_num_DoseotherControl3D',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'_run',num2str(i_sample),'.mat']);
        end
        cell_num1 = mean(cell_num1_arr,3);
        cell_num1_std = std(cell_num1_arr,[],3);
        
        cell_num2 = mean(cell_num2_arr,3);
        cell_num2_std = std(cell_num2_arr,[],3);
        cell_num3 = mean(cell_num3_arr,3);
        cell_num3_std = std(cell_num3_arr,[],3);
        cell_num4 = mean(cell_num4_arr,3);
        cell_num4_std = std(cell_num4_arr,[],3);
        
        
        
        sim_count = sim_count + 1;
        save(['Results\cell_num_Dose3D_f1920_cl3_m1_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat'],'cell_num1','cell_num1_std','tc1');
        save(['Results\cell_num_Dose3D_f1920_cl3_m2_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat'],'cell_num2','cell_num2_std','tc2');
        save(['Results\cell_num_Dose3D_f1920_cl3_m3_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat'],'cell_num3','cell_num3_std','tc3');
        save(['Results\cell_num_Dose3D_f1920_cl3_m4_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat'],'cell_num4','cell_num4_std','tc4');
        
        %clear cell_num1 cell_num1_std  cell_num1_arr
    end
end
toc