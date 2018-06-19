function[ hypoxi_level_mean, hypoxi_level_std, num_t,cell_num,o2_lev]=computehypoovertime2(fracD,numD,delT,O2_th,numsamples)
load(['Results\cellwithstacksm1m2_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'_run1.mat']);
num_t=size(o2_stack2,3); % # of time points
hypoxi_level = zeros(num_t,numsamples);
for i=1:numsamples
    load(['Results\cellwithstacksm1m2_',num2str(fracD),'nbdose',num2str(numD),'delT',num2str(delT),'_run',num2str(i),'.mat']);
    T_stack=T_stack2;
    o2_stack=o2_stack2;
    cell_num=cell_num2;
    o2_lev=o2_lev2;
    
    % figure out the amount of tissue that is hypoxic
    num_t=size(o2_stack,3); % # of time points
    
    for i_t=1:num_t
        tumor_mask = and(T_stack(:,:,i_t-1)>0,T_stack(:,:,i_t-1)~=3);% PERRINE: Tumor was was always set from tumor orginal size ((0), corrected to current step size
        hypoxi_level(i_t,i) = sum(and(o2_stack(:,:,i_t-1)<O2_th,tumor_mask))/sum(tumor_mask);
    end
end
hypoxi_level_mean=nanmean( hypoxi_level,2);
hypoxi_level_std=nanstd(hypoxi_level,[],2);
end