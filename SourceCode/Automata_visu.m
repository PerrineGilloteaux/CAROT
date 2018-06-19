function Automata_visu(cell_stack,o2_stack,o2_lev,cell_num,cell_num_std,tc,O2_th,params)
f = figure('Visible','off','Position',[360,500,200,285],'MenuBar','none');
f.Name = 'Results visualization';
b_displayNumcell    = uicontrol('Style','pushbutton',...
             'String','Plot Tumor Volume','Position',[10,220,150,25],...
             'Callback',{@b_displayNumcell_Callback});
b_export    = uicontrol('Style','pushbutton',...
             'String','Export to excel','Position',[10,180,150,25],...
             'Callback',{@b_export_Callback});
htext  = uicontrol('Style','text','String','Select Sample',...
           'Position',[10,135,40,25]);
       
hsld = uicontrol('Style', 'slider',...
        'Min',1,'Max',length(o2_stack),'Value',1,...
        'Position', [60 135 120 20],'SliderStep',[1/(length(o2_stack)-1) 1/(length(o2_stack)-1)],...
             'Callback',{@hsld_Callback}); 
b_seeO2 = uicontrol('Style','pushbutton',...
             'String','Show O2 level for s.','Position',[10,80,150,25],...
             'Callback',{@b_seeO2_Callback});
b_seeCells = uicontrol('Style','pushbutton',...
             'String','Show Tumor cell states for s.','Position',[10,40,150,25],...
             'Callback',{@b_seeCells_Callback});         
f.Visible = 'on';
s=1;
function b_displayNumcell_Callback(source,eventdata) 
% Display surf plot of the currently selected data.
figure,
sc=0.015*0.015;
shadedErrorBar(1:length(cell_num),4/3*(sc*sqrt(sc))*(sum(cell_num(:,[1,2,5]),2).*sqrt(sum(cell_num(:,[1,2,5]),2))),4/3*(sc*sqrt(sc))*sum(cell_num_std(:,[1,2,5]),2).*sqrt(sum(cell_num(:,[1,2,5]),2)));
xlabel('time in days');
ylabel('Volume in mm^3');
end


function b_export_Callback(source,eventdata) 
% Display surf plot of the currently selected data.
     header={'Normoxic_Cells',' Hypoxic_resistant_cells',...
         'Endothelial_cells', 'Healthy_Cells', 'programmed_for_mitotic_death','std_Normoxic_Cells','std_Hypoxic_resistant_cells',...
         'std_EC', 'std_Healthy_Cells', 'std_programmed_for_mitotic_death'};
     filename = 'ExportedData.xls';
     folder_name = uigetdir('.','Save results in this directory');
     ds = mat2dataset([cell_num,cell_num_std]);
     ds.Properties.VarNames =header;
    parametersfile='param.txt';
    tcfile='tc.txt';
    csvwrite([folder_name,'/',parametersfile],params);
    export(ds,'XLSfile',[folder_name,'/',filename]);
    csvwrite([folder_name,'/',tcfile],tc'); 
end

function hsld_Callback(source,eventdata) 
% Display surf plot of the currently selected data.
     s=round(get(source,'Value'));
     disp(s);
end

function b_seeO2_Callback(source,eventdata) 
% Display surf plot of the currently selected data.
     dipshow(o2_stack{s},[0.0001 0.05],jet(256));
     dipshow(o2_stack{s}<O2_th);
     figure,
     plot(o2_lev{s}*100);
     xlabel('Days')
     ylabel('Mean level of oxygen in %');
end

function b_seeCells_Callback(source,eventdata) 
% Display surf plot of the currently selected data.
     dipshow(cell_stack{s},'labels');
end
end