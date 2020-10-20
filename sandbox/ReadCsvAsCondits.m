%% ReadCsv  -  alpha / beta
%   take path and read in each csv
%   return exper and condits
%
%       exper
%           conditions
%           conditWellMap
%           groupConditMap
%           conditIndexMap
%           frames
%
%       condits
%           name
%           wells
%               name
%               path (csv)
%               raw
%               cellCount
%               cells
%                   <COLS>
%



function [exper, condits] = ReadCsvAsCondits(experPath)
    global CONST
    [COL_LAYOUT,COL_COUNT,COLS] = AeplUtil.GetColLayout(CONST.COL_LAYOUT_VER); 
     
    
    csvDir = fullfile(experPath, CONST.CSV_DIR);
    
   
    exper = struct();
    disp(experPath)
    disp(CONST.PLATE_MAP_SUF)
    plateMapFile = AeplUtil.FindFile(experPath,CONST.PLATE_MAP_SUF);
    %if plateMapFile
        % could all be done in one line but it was really long
            [a,b,c] = ReadPlateMap(plateMapFile);
            exper.conditions = a;
            exper.conditWellMap = b;
            exper.groupConditMap = c;
        
    %end
    
    
    clear condits       % is this necessary?
    % use conditWellMap and conditIndexMap to allow condition names to not
    % be valid matlab identifiers
    %       I do this by returning from ReadPlateMap, an
    exper.conditIndexMap = containers.Map();
    condits(length(exper.conditions)) = struct();
    for c = 1:length(exper.conditions)
        condits(c).name = exper.conditions{c};  
        exper.conditIndexMap(exper.conditions{c}) = c; 

        disp(condits(c).name)
        disp(exper.conditWellMap(exper.conditions{c}))

        conditWells = exper.conditWellMap(exper.conditions{c});
        condits(c).wells(size(conditWells)) = struct();

        for w = 1:length(conditWells)
            condits(c).wells(w).name = conditWells{w};
            %condits(c).wells(w).path = strcat(exper.folder,'Xls/', exper.name,'_',condits(c).wells(w).name{1}, '.csv');
            condits(c).wells(w).path = AeplUtil.LookForFile(csvDir,condits(c).wells(w).name);
            %if exist(condits(c).wells(w).path, 'file')
            if condits(c).wells(w).path

                
                    
                condits(c).wells(w).raw = readtable(condits(c).wells(w).path);
                    
                condits(c).wells(w).cellCount = width(condits(c).wells(w).raw)/COL_COUNT;

                condits(c).wells(w).raw = condits(c).wells(w).raw.Variables;

                %% Read in csv file
%                 condits(c).wells(w).raw = readtable(condits(c).wells(w).path);
%                 condits(c).wells(w).raw = condits(c).wells(w).raw.Variables;
%                 condits(c).wells(w).cellCount = width(condits(c).wells(w).raw)/COL_COUNT;
                
                %% Intialize struct array for cells
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
                fprintf(1,'missing csv file for well %s.\n', condits(c).wells(w).name)
            end
        end
    end
    
    exper.frames = 0;
    % let's make this frameCount
    while (exper.frames == 0)
        try 
            exper.frames = length(condits(1).wells(1).cells(1).xcoords);
        catch e
        end
    end

end