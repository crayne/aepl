%% RunSummarizeData  -  alpha
%   wrapper function for heavy-lifters to read in data and make
%   csv files with displacement data by condition and by well
%
%   requires plate map, but doesn't require groups
%


function RunSummarizeData(experPath,exper,condits)
    
    %[exper,condits] = ReadCsvAsCondits(experPath);
    

    global CONST
    procDir = fullfile(experPath, CONST.PROCESSED_DIR);
    if ~exist(procDir,'dir')
        mkdir(procDir)
    end  
    
    exper.t_int = (1/6):(1/6):(exper.frames+10/6);

    t_int = num2cell(exper.t_int(1:exper.frames-1)');
    
    
    
    
    conditAvgC = cell(0);
    conditMediC = cell(0);
    
    wellAvgC = cell(0);
    wellMediC = cell(0);
    
    caCol = 1;
    cmCol = 1;
    
    waCol = 1;
    wmCol = 1;
    
        
    conditAvgC{1,caCol} = 'Time (hours)';
    conditAvgC(2:exper.frames,caCol) = t_int;
    caCol = caCol + 1;
    
    conditMediC{1,cmCol} = 'Time (hours)';
    conditMediC(2:exper.frames,cmCol) = t_int;
    cmCol = cmCol + 1;
    
    
    wellAvgC{1,waCol} = 'Time (hours)';
    wellAvgC(2:exper.frames,waCol) = t_int;
    wCaol = waCol + 1;
    
    wellMediC{1,wmCol} = 'Time (hours)';
    wellMediC(2:exper.frames,wmCol) = t_int;
    wmCol = wmCol + 1;
    
    
    for condit = condits
        
        %if ~isfield(condit,mat) 
        condit  = AeplUtil.MakeConditMat(condit);
        %end
       
        conditAvg = nanmean(condit.mat,2);
        conditAvgC{1,caCol} = condit.name;
        conditAvgC(2:length(conditAvg)+1, caCol) = num2cell(conditAvg);
        caCol = caCol + 1;
        
        
        conditMedi = nanmedian(condit.mat,2);
        conditMediC{1,cmCol} = condit.name;
        conditMediC(2:length(conditMedi)+1, cmCol) = num2cell(conditMedi);
        cmCol = cmCol + 1;
        
%         for w = 1:length(condit.wells)
%             wellAvg = nanmean(condit.wells(w).mat,2);
%             wellAvgC{1,waCol} = condit.wells(w).name;
%             if ~isempty(wellAvg)
%                 wellAvgC(2:length(wellAvg)+1, waCol) = num2cell(wellAvg);
%             end
%             waCol = waCol + 1;
%             
%             wellMedi = nanmean(condit.wells(w).mat,2);
%             wellMediC{1,wmCol} = condit.wells(w).name;
%             if ~isempty(wellMedi)
%                 wellMediC(2:length(wellMedi)+1, wmCol) = num2cell(wellMedi);
%             end
%             wmCol = wmCol + 1;
%         end
    end
    
    
 
    conditT = cell2table(conditAvgC);
    conditFile = fullfile(procDir,[CONST.AVG_PREF, CONST.CONDIT_DATA_SUF]);
    writetable(conditT, conditFile, 'WriteVariableNames',0)
    
    conditT = cell2table(conditMediC);
    conditFile = fullfile(procDir,[CONST.MEDI_PREF, CONST.CONDIT_DATA_SUF]);
    writetable(conditT, conditFile, 'WriteVariableNames',0)
    

    wellT = cell2table(wellAvgC);
    wellFile = fullfile(procDir,[CONST.AVG_PREF, CONST.WELL_DATA_SUF]);
    writetable(wellT, wellFile, 'WriteVariableNames',0)
    
    wellT = cell2table(wellMediC);
    wellFile = fullfile(procDir,[CONST.MEDI_PREF, CONST.WELL_DATA_SUF]);
    writetable(wellT, wellFile, 'WriteVariableNames',0)



end