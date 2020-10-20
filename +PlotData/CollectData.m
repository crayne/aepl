function [data,well] = CollectData(root)

dlist = dir(root);

data = cell(length(dlist),1);
well = cell(length(dlist),4);
for i = 1:length(dlist)
    name = fullfile(dlist(i).folder,dlist(i).name);
    [~,f,~] = fileparts(dlist(i).name);
    well{i,1} = f(1);
    well{i,2} = str2double(f(2:end));
    well{i,3} = double(f(1)-64);
    CellLine = regexp(dlist(i).folder, filesep, 'split');
    
    well{i,4} = CellLine{3};
    well{i,5} = CellLine{4};
    data{i} = readtable(name);
end

end 

