%
% Simulation of the Prostate tumor irradiation
% Model assumes the following:
% radiation protocol defined by: # of exposure (numD), dose per fraction
% (fracD), time between exposure (delT), % of radioresistant cells RTres
% leak_death: % (>1) each time there is a death hit in blood vessel, that vessel becomes leak_death time more leaky
% extra_days: % number of days to follow after radiation (default: 30)
%
%
% Tumor Cell Properties
% Tumor cells have blood vessel intercalated with O2 leakiness parameter
% Tumor cells have an oxygen requirement, if low no proliferation for tumor
% type 1 (need O2). Tumor type 2 are initially lower in frequency but do
% not need O2 and are radioresistant by a factor resist_factor.
%
% Radiation effect
% 1. Death is a function of O2: dose equ = (0.6*perc_O2 + 0.4)*fracD
% 2. Radiation makes blood vessel leakier for O2 - D dependence unknown
% 3. Radiation kills blood vessel cells - Use same as rule 1 for dose
% dependence - No blood cell means no more O2
%
%
% Sylvain V. Costes
% November 2014
% Berkeley Lab
% Version 4: Important modification from version 3.
% Assume cell death occurs only when a cell labeled to die enter mitosis (mitotic catastrophe)
% Version 5: (no perf leak factor has to be set to 1 instead of 1.3-> no difference with version 4)
% NO O2 effect: HRF set to FracD as an option O2 effect false.
% No vessel Death: vessel death turned off is option EC_death_option set to false
% Version 6: Done by Perrine
% Version 7: Added # of Cell layers that can divided as a parameter for division
%
% Version 10: adding a waitfor variable not to start the irradiation from
% the beginning
% Version 11: Modifying the rules for hypoxia: can no divide as others +
% change the connexity (to avoid squaure tumors) with elliptic struring
% element by using dilation instead of bdilation (and consequently multiplied by 2 the size, which corresponds to
% a rdius while it was the diameter in bdilation
%changed the size of the tumor ratio to image to avoid saturation
function [cell_stack,o2_stack,cell_num,o2_lev] = Hypoxia_modelvpp(field_size,protocol,RTres,O2_th,leak_death,extra_days,O2_option,EC_death_option,cell_layer,parralel,alpha,beta,alv,bev,blood_density)

%% Create tumor tissue
if ~exist('extra_days')
    extra_days = 30;
end
if ~exist('cell_layer')
    cell_layer = 1;
end
if ~exist('blood_density')
blood_density = 0.037963;
end
vessel_orientation = 'z';
%blood_density = 0.00922; % proportion of cell being blood - to reach 1.3% average O2 level in tumor
blood_density = 0.037963; % proportion of cell being blood - to get 1% hypoxic region only
dif_cst = 2.4; % diffusion constant used for oxygen. Obtained from another simulation (O2_diffusion.m)
Pscale = 1.18; % Value at 15 um away to scale on.
Pcoef=[1.31e6,-3.12e4,6.58e2,0.5]; % polynomial used to describe HRF
cell_array = ones(field_size,field_size); % the tissue - 1 represents a normal cell
tumor_mask = newim(field_size,field_size);
tumor_mask(round(end/2),round(end/2)) = 1;
tumor_mask = dilation(tumor_mask,field_size/6);
cell_array = double(dip_image(cell_array)*(tumor_mask>0));
num_cell = sum(tumor_mask>0);

% Create blood vessell for various traversal geometry
rand_array = rand(field_size, field_size);
switch vessel_orientation
    case 'x'
        cell_array(field_size/20:field_size/5:end,:) = 3; % blood vessel are type 3.
    case 'y'
        cell_array(:,field_size/20:field_size/5:end) = 3; % blood vessel are type 3.
    case 'z'
        cell_array(rand_array<blood_density)= 3; % blood vessel are type 3.
end
% Death based on clonogenic data of Vincent from PC3 tumor cells
if ~exist('alpha')
    alpha = 0.0441;
end
if ~exist('beta')
    beta = 0.0898;
end

% Death of vessels based on Riquier 2013 Rad Ther and Onc
if ~exist('alv')
   alv = 0.19;
end
if ~exist('bev')
   bev = 0.039;
end
% O2 levels below which type 1 cells cannot divide
if ~exist('O2_th')
    O2_th = 0.002;
end
lengthprotocol=length(protocol);
% Growth arrest duration (days) after radiation - PIECEWISE cubic interpolation -
% Data to be updated by Francois team
d_ar=[0 1 2 4 8 15]; %dose for growth arrest measurement
gar = [0 0 24 36 48 65]/24;

 fracD=protocol(1); % coud be 9, then gar_days will be 0
 gar_days = interp1(d_ar,gar,fracD,'pchip'); %for growth arrest to heppen from first day

    % Number of days simulated
sim_days = lengthprotocol*1.2;

Tcycle = 1; % cell cycle time in day
timer_day = 0; % Keep track of time in days
timer_rad = 0; % PERRINE: WARNING: start at delT-1 instead of 0: ortherwise first week would have a different protocol; % keep track of interval between radiation exposure
dose_cnt = 0;
num_step = fix(sim_days/Tcycle) + extra_days; % step at each cell cycle time (easier this way) - Add extra days to see if tumor grows back
nprotocol=zeros(num_step,1);
nprotocol(1:length(protocol))=protocol;
protocol=nprotocol;
if parralel==0
    cell_stack = newim(field_size,field_size,num_step); % keep snapshot of cell types per time step
    o2_stack = newim(field_size,field_size,num_step); % keep snapshot of O2 level per time step
else
    cell_stack=[]; %% to speed up computation since can not be saved
     o2_stack =[];
end

cell_num = zeros(num_step,5); % keeps track of cell number as a function of time. All dead cells are in index 5
o2_lev = zeros(num_step,1); % keeps track of oxygen levels as a function of time
% gar_array keeps track of growth arrested cell.
% If value equal gar_days, cell are not arrested. If lower, increment by
% time and do not divide

% Populate with hypoxic resistant tumor cells (ID=2)
vessel_img = 1.0*dip_image(cell_array==3); % this initializes position of blood vessel all to 1. As they get damaged, this value increases.
perc_O2 = min(max(gaussf(vessel_img*Pscale,dif_cst),0.001),0.05);
rt_mask = and(tumor_mask>0,perc_O2<O2_th);
num_rt = round(RTres*num_cell); % number of radioresistant cells
rt_ind = find(rt_mask>0); % index showing which pixel is in hypoxia
% try % if it does not work, will cluster the hypoxic cells into one region
%     rt_ind = randswap(rt_ind); % randomize order
% end
% rt_mask(rt_ind(num_rt+1:end))=0; % get rid of all extra cells
%cell_array(double(rt_mask)>0) = 2;
cell_array(double(rt_ind)) = 2;
% Populate rest of tissue with healthy cells, with 20% occupancy (Leave
% some gaps)
cell_array(and(rand_array<0.2,cell_array==0)) = 4;


%% Start loop
for i_t = 1:num_step % at each step (corresponding to a 24h cell division)
    cell_num(i_t,:) = [sum(cell_array(:)==1),sum(cell_array(:)==2),sum(cell_array(:)==3),sum(cell_array(:)==4),sum(cell_array(:)>4)];
    cell_img = dip_image(cell_array);
cell_img(cell_img>4) = cell_img(cell_img>4) + 1; % Increment all dying cells by one day.
    %    cell_img(cell_img==5+T_clear) = 0; % Remove all dead cells that have been in tissue for more than T_clear days.
    perc_O2 = min(max(gaussf(vessel_img*Pscale,dif_cst),0.001),0.05);
 
    if  parralel==0
        cell_stack(:,:,i_t-1) = cell_img*(cell_img<5) + 5*(cell_img>=5); % save cells, all dying cells set to 5
        o2_stack(:,:,i_t-1) = perc_O2;
    end
    timer_day = timer_day + Tcycle;
    timer_rad = timer_rad + Tcycle;
    
    try
        o2_lev(i_t) = mean(perc_O2(and(tumor_mask,cell_img<3))); % O2 levels in tumor only
    catch
        o2_lev(i_t) = 0;
        break;
    end
       if timer_rad>gar_days
       gar_days=0;
        % refill array with adjacent dividing cells. Radioresistant type 2
        % cells grow better in low O2
        % Cells set to die can also divide and will be removed if selected
        % to divide
        fill_index = and(cell_img==0,dilation(or(or(cell_img==1,cell_img==2),cell_img>4),cell_layer*2)); % space to be filled - empty pixels cell_layer away from a cell that can divide
       % div_index = and(or(or(and(cell_img==1,perc_O2>O2_th),cell_img==2),cell_img>4),bdilation(fill_index,cell_layer)); % potential cells that can divide
       % ALl cells can now divude eve in hypoxia 
       div_index = and(or(or(cell_img==1,cell_img==2),cell_img>4),dilation(fill_index,cell_layer*2)); % potential cells that can divide
        roi = or(div_index,fill_index); % Region of interest including dividing cells and open space for them
        % fill dividing cells layer by layer (otherwise, outerlayer endup
        % duplicating more than once for cell_layer>1
        left_img = fill_index; % set initially to have all fill_index available. Reduce as it gets filled by new divided cell for each layer dividing
%        visual_layer = newim(size(cell_img)); % for debug to visualize
%        growth
        for i_layer = 1:cell_layer
            layer_index = and(div_index,dilation(fill_index,i_layer*2)); % cells to consider for division in this layer.
            max_div = sum(layer_index>0); % # of new cells should be limited to available space and number of cells that can divide
            div_index = xor(div_index,layer_index); % Remove the cell that were just considered for division in this layer.
            % use larger dilation as as we move away from the first layer,
            % deeper layers cannot reach left_img if this multiplicative
            % factor 2 is not there...
            spread_img = left_img*dip_growregions(uint8(cell_img*layer_index),[],roi,2,round(i_layer*2),'low_first'); % Let dividing cells expand using watershed of radius equal to i_layer*2 then only keep the one in the spread_img boolean
            spread_ind = find(spread_img>0);
            num_fill = size(spread_ind(:),1);
            % only keep max_div fill gap randomly to make sure division only
            % multiplies by 2
            spread_ind = spread_ind(randperm(num_fill));
            spread_img(spread_ind(max_div+1:num_fill)) = 0; % make sure only max_div cell dividing. Set the rest to 0
%            left_img = xor(left_img,spread_img>0);
%           MAke sure that all labels 2are set back to type 1 (to avoid type 2
%           division)
             spread_img(spread_img==2)=1;
             
            cell_img(spread_img>0) = spread_img(spread_img>0);
            left_img = and(fill_index,cell_img==0); % remove from fill_index pixels now filled with cell
            % This is for debug to visualize growth
%             visual_layer(spread_img>0) = i_layer;
%             dipshow(visual_layer(100:300,230:280),'labels');
%             diptruesize(600)
        end
        % now cells can divide even in hypoxic regions
        %cell_img(and(cell_img==1,perc_O2<O2_th))=0; % All type 1 cells cannot be in hypoxic region and if they were created during division, they are removed
   
        cell_img(and(cell_img>4,roi)) = 0; %Any cells committed to die that have divided is killed.
        %cell_img = cell_img + spread_img*fill_index; %replaces all the cells that have been killed with the new divided cells
    end
    cell_array = double(cell_img);
   if protocol(i_t)~=0 % personalized protocol indicate a radiation
        fracD=protocol(i_t);% Time to irradiate, except on day 6 and 7 (weekend)
       mydays='LMmJVSD';
       disp(mydays(mod(i_t,7)));
       gar_dayst = interp1(d_ar,gar,fracD,'pchip'); 
       if gar_dayst>gar_days-timer_rad
           gar_days=gar_dayst
           timer_rad = 0; 
       end
        if O2_option
            hrf_fit = 1+1./polyval(Pcoef,double(perc_O2));
            Deq = fracD./hrf_fit; % Equivalent dose based on O2 level
            % Probability of death for each fractionated dose
            Pdeath = 1-exp(-alpha.*Deq-beta.*Deq.^2); % This is the probability to kill a cell
        else
            Pdeath = 1-exp(-alpha*fracD-beta*fracD^2);
        end
        dose_cnt = dose_cnt + 1;
        
        Pdeath_v = 1-exp(-alv*fracD-bev*fracD^2); % This is the probability to hit a vessel and increase vessel leakiness
        Pdeath_vessel = arrayfun(@vessel_death,fracD); % This is the probability to remove a vessel from very high dose
        rand_array = rand(field_size,field_size);
        % kill cells with frequency beta - This is done by indexing them to
        % 5. Then that number is incremented at each time step until it
        % reach T_clear. This is when cell is removed
        cell_array(and(and(Pdeath>rand_array,cell_array<3),cell_array>0)) = 5;
       
        vessel_img(dip_image(and(rand_array<Pdeath_v,cell_array==3))) = vessel_img(dip_image(and(Pdeath_v>rand_array,cell_array==3)))*leak_death;
        % To be turned off if no vessel death
        if EC_death_option
            vessel_img(dip_image(and(rand_array<Pdeath_vessel,cell_array==3))) = 0;
            cell_array(and(rand_array<Pdeath_vessel,cell_array==3))= 0;
        end
         % set the cell with hypoxia indicator
       
   end
        % Let cell divide only if tissue is not growth arrested

     cell_array(and(or(cell_array==1,cell_array==2),perc_O2<O2_th))=2;
end


