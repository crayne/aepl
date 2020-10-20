%% RunPlotData  -  beta
%   wrapper function for calling the heavy-lifters that make plots
%
%   assumes plate map has groups
%


function RunPlotData3(experPath)
    

    global CONST
    

    [exper,condits] = ReadCsvAsCondits(experPath);
    

    procDir = fullfile(experPath, CONST.PROCESSED_DIR);
    if ~exist(procDir,'dir')
        mkdir(procDir)
    end    

    exper.t_int = (1/6):(1/6):(exper.frames+10/6);
    
    exper.ylimit = 20;
    
    groups = exper.groupWellMap.keys();
   
    condits2(length(condits)) = struct('name','','wells',[],'mat',[]);

    %temp = exper.conditIndexMap('DMSO%');
    temp = exper.conditIndexMap('RH30, DMSO .1%');
%              try
    dmso = AeplUtil.MakeConditMat(condits(temp));
    
    %condits2(cond)= temp;
    
%     if plotOptions.laeoutSet.make 
%         laeout = CONST.SET_LAEOUT;
%         count = 0;
%     % elseif plotOptions.laeoutByOneGroup
%         % need to add
%     end
    
    for group = groups
        if strcmp(group{1},'DMSO')
            continue
        end
%         if plotOptions.laeoutSet.make 
%             count = 0;
        gCondits = exper.groupWellMap(group{1});
        
%         if plotOptions.laeoutByGroup.make
             plotCount = length(gCondits);
             laeout = AeplUtil.MakeLaeout(plotCount);
%         end
        
        figure
        plotNum = 1;
        for conditNamee = gCondits
            conditName = conditNamee{1};
            cond = exper.conditIndexMap(conditName);
%              try
               temp = AeplUtil.MakeConditMat(condits(cond));
               condits2(cond)= temp;
%             catch e
%                 disp(conditName)
%                 disp(e.getReport)
%             end
            subplot(laeout(1),laeout(2), plotNum)
            title(conditName)
            
            PlotData.MakeConditSubplot2(condits2(cond).mat,exper,dmso)
            plotNum = plotNum + 1;
            
%             if plotOptions.laeoutSet.make 
%                 if plotNum > laeout(1)*laeout(2)
%                     h = gcf;
% 
%                     set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')
% 
%                     fname = [group{1},'_',num2str(count),CONST.PLOT_SUF];
%                     print(h,fullfile(procDir,fname),'-dpdf')
%                     close(h)
%                   
%                     figure
%                     count = count+1;
%                     plotNum = 1;
%                 end
%             end
        end
%         if plotOptions.laeoutSet.make 
           if plotNum > 1
                h = gcf;

                set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')

                fname = [group{1},CONST.PLOT_SUF];

                print(h,fullfile(procDir,fname),'-dpdf')
                close(h)
           else
               close(gcf)
           end
%         else
%             h = gcf;
% 
%             set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')
% 
%             fname = [group{1},CONST.PLOT_SUF];
% 
%             print(h,fullfile(procDir,fname),'-dpdf')
%             close(h)
% 
%         end
        
    end
    %condits = condits2;
end
    
