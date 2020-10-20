%%RunPlotData
clear

root = 'G:\MSKCC_Baylies-Lab_Rhabdo-Data\*\*\*map.csv';
dlist = dir(root);

%%
bremove = false(size(dlist));
for i = 1:length(dlist)
    name = dlist(i).name;
    if name(1)=='.'
        bremove(i) = 1;
        
    end
end
dlist(bremove) = [];

for i = 1:length(dlist)
    name = fullfile(dlist(i).folder,dlist(i).name);
    
    [~,~,D] = xlsread(name);
    meta{i} = D;
    %%
end

%%
root = 'G:\MSKCC_Baylies-Lab_Rhabdo-Data\*\*\Czi\Csv\*.csv';

%%

[rawdata,well] = PlotData.CollectData(root);

profile on 
Data = PlotData.MakeFeatures(rawdata);
profile viewer

featureNames = {'Mitosis','Death','Death %','Population','Migration';
                '(# Cells)','(# Cells)','(% Cells)','(# Cells)','(cm/min)'}';
%%
close all

PlotData.ExportData(Data,well,featureNames);


PlotData.MakeFigures(Data,well,featureNames)













