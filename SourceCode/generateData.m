function [cell_stack,o2_stack,o2_lev,cell_num,cell_num_std,tc,computetime] = generateData(field_size,fracD,numD,delT,RTres,O2_th,leak_death,extra_days,O2_option,EC_death_option,n_cell_layer,parallel,num_sample,alpha,beta,alv,bev,blood_density )
tic

tc = nan(num_sample,1);



%% Data for growth control only

fprintf('%d Gy, %d fractions, %d days between exposure\n',fracD,numD,delT);
cell_num_arr=[];
cell_num=[];
leaking_factor=leak_death;

    
%% actual simulation for numcell only (generate stack where needed with parralel option set to 0 and a classical for loop
% (one patient is enough for stack by the way)






parfor i_sample = 1:num_sample; % number of repeat for average
    
    % Tstack not used here for parfor: use a classical for loop to
    % get T_stack1
    [cell_stack{i_sample},o2_stack{i_sample},cell_num,o2_lev{i_sample}] = Hypoxia_model(field_size,fracD,numD,delT,RTres,O2_th,leaking_factor,extra_days,O2_option,EC_death_option,n_cell_layer,parallel,alpha,beta,alv,bev,blood_density);
    
    tc(i_sample) = sum(cell_num(end,1:2))==0;
    cell_num_arr(:,:,i_sample) = cell_num;
    
    
    % forbidden in parfor save(['Results\cell_num_DoseotherControl3D',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'_run',num2str(i_sample),'.mat']);
end
cell_num = mean(cell_num_arr,3);
cell_num_std = std(cell_num_arr,[],3);

%save(['Results\cell_num_Dose3D_f1920_cl3_m1_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'test2.mat'],'cell_num','cell_num_std','tc','T_stack','o2_stack','o2_lev');



computetime= toc;

