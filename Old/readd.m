function [exper,condits] = readd(mypath)

COL_COUNT = 5;

COL_LAYOUT.xcoords = 1;
COL_LAYOUT.ycoords = 2;
COL_LAYOUT.area = 3;
COL_LAYOUT.covarience = 4;
COL_LAYOUT.distance = 5;
COLS = (fieldnames(COL_LAYOUT));

PLATE_MAP_SUF = '_plate-map.csv';


exper = struct();
exper.folder = mypath;

[~,exper.name] = fileparts(exper.folder(1:end-1));
exper.plateMapFile = strcat(exper.folder,exper.name,PLATE_MAP_SUF);

%[exper.conditions,exper.conditWellMap, exper.groupWellMap] = plateMap2condits(exper.plateMapFile);
[exper.conditions,exper.conditWellMap, exper.groupWellMap] = readPlateMap(exper.plateMapFile);


%dlist2 = dir(fullfile(exper.folder,'Xls','*.csv'));
clear condits
exper.conditIndexMap = containers.Map();
condits(size(exper.conditions,2)) = struct();
for c = 1:size(exper.conditions,2)
    condits(c).name = exper.conditions{c};  
    exper.conditIndexMap(exper.conditions{c}) = c; 
    
    disp(condits(c).name)
    disp(exper.conditWellMap(exper.conditions{c}))
    
    conditWells = exper.conditWellMap(exper.conditions{c});
    condits(c).wells(size(conditWells)) = struct();

    for w = 1:size(conditWells,2)
        condits(c).wells(w).name = conditWells(w);
        condits(c).wells(w).path = strcat(exper.folder,'Xls/', exper.name,'_',condits(c).wells(w).name{1}, '.csv');
        if exist(condits(c).wells(w).path, 'file')
            
            condits(c).wells(w).raw = readtable(condits(c).wells(w).path);
            
            condits(c).wells(w).cellCount = width(condits(c).wells(w).raw)/COL_COUNT;

            condits(c).wells(w).raw = condits(c).wells(w).raw.Variables;

            condits(c).wells(w).cells(condits(c).wells(w).cellCount) = struct();
           
            for j = 1:condits(c).wells(w).cellCount
                startCol = (j-1)*5;
                condits(c).wells(w).cells(j).id = condits(c).wells(w).raw{1,startCol+1};
                
                for k = 1:COL_COUNT
                    condits(c).wells(w).cells(j).(COLS{k}) = str2double(condits(c).wells(w).raw(5:end,startCol+COL_LAYOUT.(COLS{k})));
                    if isnan(condits(c).wells(w).cells(j).(COLS{k})(end))
                        condits(c).wells(w).cells(j).(COLS{k})(end) = [];
                    end
                end
               
            end
        else 
            disp(strcat('missing csv file', condits(c).wells(w).name))
        end
    end
end

end