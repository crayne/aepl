function wellAndConditAvgs(exper, condits)
    SUB_DIR = 'summaryCsv';
    CONDIT_SUF = '_Avg-veloc-conditions.csv';
    WELL_SUF = '_Avg-veloc-wells.csv';
    
    



    conditC = cell(0);
    wellC = cell(0);
    
    cCol = 1;
    wCol = 1;
    
    
    t_int = num2cell(exper.t_int(1:exper.frames-1)');
    
    conditC{1,cCol} = 'Time (hours)';
    conditC(2:exper.frames,cCol) = t_int;
    cCol = cCol + 1;
    
    wellC{1,wCol} = 'Time (hours)';
    wellC(2:exper.frames,wCol) = t_int;
    wCol = wCol + 1;
    
    for condit = condits
        
        if ~isfield(condit,mat) 
            condit.mat = 
        conditAvg = nanmean(condit.mat,2);
        conditC{1,cCol} = condit.name;
        conditC(2:length(conditAvg)+1, cCol) = num2cell(conditAvg);
        cCol = cCol + 1;
        
        for w = 1:length(condit.wells)
            wellAvg = nanmean(condit.wells(w).mat,2);
            wellC{1,wCol} = condit.wells(w).name;
            if ~isempty(wellAvg)
                wellC(2:length(conditAvg)+1, wCol) = num2cell(wellAvg);
            end
            wCol = wCol + 1;
        end
    end
    
    
    %[P,F,~] = fileparts(Name);
    myPath = fullfile(exper.folder, SUB_DIR);
    if ~exist(myPath,'dir')
        mkdir(myPath)
    end
    
    conditT = cell2table(conditC);
    conditFile = fullfile(myPath,[exper.name CONDIT_SUF]);
    writetable(conditT, conditFile, 'WriteVariableNames',0)
    
    wellT = cell2table(wellC);
    wellFile = fullfile(myPath,[exper.name WELL_SUF]);
    writetable(wellT, wellFile, 'WriteVariableNames',0)



% T = cell2table(Mat);
% writetable(T,XlsName,'WriteVariableNames',0);