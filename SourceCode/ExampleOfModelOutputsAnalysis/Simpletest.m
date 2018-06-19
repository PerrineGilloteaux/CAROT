
field_size = 400; %tumor dilation is fieldsize/4 not by 2 to avoid space saturation


extra_days = 15; % number of days to track after IR (important to see if cell goes to 0)
RTres=0.01; % At beginning, difonly 1% of tumor cells are hypoxic.
O2_th=0.002; % 0.2% Oxygen level defines hypoxia below it.

%sim_results = zeros(11*4,6); % For each dose and fraction, report fracD, numD, TotalD, Tumor control, TC error, delT
% 4 doses, 11 fraction regimens


sim_count = 1;
T_stack1=[];
o2_stack1=[];


o2_lev1=[];
indexD=0;

n_cell_layer=3;
fracD = 2;
numD=1;
delT=1;
   
        
fprintf('%d Gy, %d fractions, %d days between exposure\n',fracD,numD,delT);
cell_num1_arr=[];


[T_stack1,o2_stack1,cell_num1,o2_lev1] = Hypoxia_model(field_size,fracD,numD,delT,RTres,O2_th,1.3,extra_days,1,1,n_cell_layer,0);

cell_num1_arr(:,:) = cell_num1;


cell_num1 = mean(cell_num1_arr,3);
cell_num1_std = std(cell_num1_arr,[],3);

TumorMask=(T_stack1==1)+(T_stack1==2)+(T_stack1==5)
T_stackmoreinfo=T_stack1+5*and(TumorMask>0,o2_stack1<O2_th)
for(t=0:length(cell_num1)-1)
    nbhypoxy(t+1)=sum(sum(and(TumorMask(:,:,t)>0,o2_stack1(:,:,t)<O2_th)));
end
