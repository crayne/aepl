%% RunPlotData  -  beta
%   wrapper function for calling the heavy-lifters that make plots
%
%   assumes plate map has groups
%


function RunPlotData4(experPath,exper,condits)
    

    global CONST
    

    %[exper,condits] = ReadCsvAsCondits(experPath);
    

    procDir = fullfile(experPath, CONST.PROCESSED_DIR);
    if ~exist(procDir,'dir')
        mkdir(procDir)
    end    

    exper.t_int = (1/6):(1/6):(exper.frames+10/6);
    
    exper.ylimit = 20;
    
    groups = exper.groupConditMap.keys();
   
    %condits2(length(condits)) = struct('name','','wells',[],'mat',[]);

    %temp = exper.conditIndexMap('DMSO%');
    temp = exper.conditIndexMap('RH30, DMSO .1%');
%              try
    dmso = AeplUtil.MakeConditMat(condits(temp));
    dmso.mat(dmso.mat>50) = nan;
    dmsoMatRowMeds = nanmedian(dmso.mat,2);
    t_int_temp2 = permute(exper.t_int(1:size(dmsoMatRowMeds,1)),[2 1]);
    %condits2(cond)= temp;
    
%     if plotOptions.laeoutSet.make 
%         laeout = CONST.SET_LAEOUT;
%         count = 0;
%     % elseif plotOptions.laeoutByOneGroup
%         % need to add
%     end
    %laeout = [2,3];
    laeout = [1,2];
    plotNum = 1;
    figNum = 1;
    
    figure
    black = [0,0,0];
    %white = [255,255,255];
    white = [.8,.8,.8];
    red = [1, 0, 0];

    %len = length(groups) - 1;
    %colors_p = [linspace(red(1),white(1),len)', linspace(red(2),white(2),len)', linspace(red(3),white(3),len)'];
    for group = groups
        if strcmp(group{1},'Control')
            continue
        end
%         if plotOptions.laeoutSet.make 
%             count = 0;
        gCondits = exper.groupConditMap(group{1});
        
        len = length(gCondits);
        colors_p = [linspace(black(1),white(1),len)', linspace(black(2),white(2),len)', linspace(black(3),white(3),len)'];
    
        
%         if plotOptions.laeoutByGroup.make
             %plotCount = length(gCondits);
             %laeout = AeplUtil.MakeLaeout(plotCount);
%         end
        
        
        %figure
        %if 
        hold on 

        subplot(laeout(1),laeout(2),plotNum);

        scatter(t_int_temp2,dmsoMatRowMeds, 10,'filled','MarkerFaceColor',red) ;
        %plotNum = 1;
        %coun = 1;
        %legendNames = cell(length(gCondits));
        %legendNames(coun) = {'DMSO'};
        %coun = coun + 1;
        subCount = 1;
        for conditNamee = gCondits
            conditName = conditNamee{1};
            %legendNames(coun) = conditNamee;
            %coun+=1
            cond = exper.conditIndexMap(conditName);
%              try
           temp = AeplUtil.MakeConditMat(condits(cond));
           %condits2(cond)= temp;
           mat = temp.mat;
           mat(mat>50)=nan;
    %matRowMeans = nanmean(mat,2);
            matRowMeds =  nanmedian(mat,2);
            t_int_temp = permute(exper.t_int(1:size(matRowMeds,1)),[2 1]);
            hold on
            scatter(t_int_temp,matRowMeds, 10,'filled','MarkerFaceColor',colors_p(subCount,:)) ;
            subCount = subCount + 1;

%             catch e
%                 disp(conditName)
%                 disp(e.getReport)
%             end
            %subplot(laeout(1),laeout(2), plotNum)
            %title(conditName)
            
            %PlotData.MakeConditSubplot2(condits2(cond).mat,exper,dmso)
            %plotNum = plotNum + 1;
            
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
        
        ylim([0 exper.ylimit])
    
        xlabel('time (hours)')              % x-axis label
        ylabel('displacement (pixels)')
        temp = [ { 'DMSO' }, gCondits ];
 
        legend(temp)
        title(group)
        plotNum = plotNum + 1;
%         if plotOptions.laeoutSet.make 
%            if plotNum > 1
%                 h = gcf;
% 
%                 set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')
% 
%                 fname = [group{1},CONST.PLOT_SUF];
% 
%                 print(h,fullfile(procDir,fname),'-dpdf')
%                 close(h)
%            else
%                close(gcf)
%            end
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
        if plotNum > laeout(1)*laeout(2)
            h = gcf;

            set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')

            fname = [num2str(figNum),CONST.PLOT_SUF];
            print(h,fullfile(procDir,fname),'-dpdf')
            close(h)

            figure
            figNum = figNum+1;
            plotNum = 1;
        end


        
        
    end
    
    if plotNum > 1
        h = gcf;

        set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')

        fname = [num2str(figNum),CONST.PLOT_SUF];

        print(h,fullfile(procDir,fname),'-dpdf')
        close(h)
    else
       close(gcf)
    end
    
    %condits = condits2;
end
    
