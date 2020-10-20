%% RunPlotData  -  beta
%   wrapper function for calling the heavy-lifters that make plots
%
%   assumes plate map has groups
%


function RunPlotData5(experPath,exper,condits)
%function RunPlotData5(experPath)
    %[exper,condits] = ReadCsvAsCondits(experPath);

    global CONST

    procDir = fullfile(experPath, CONST.PROCESSED_DIR);
    if ~exist(procDir,'dir')
        mkdir(procDir)
    end    

    exper.t_int_long = (1/6):(1/6):(exper.frames+10/6);
    
    exper.ylimit = 20;
    
    groups = exper.groupConditMap.keys();
   
    
    %temp = exper.conditIndexMap(CONST.CONTROL);
    
%     tempK = AeplUtil.FindKey(exper.conditIndexMap,CONST.CONTROL);
%     if tempK{1}
%         temp = exper.conditIndexMap(tempK{1});
%     else
%         %error('no control')
%     end

    temp = AeplUtil.useKeyPattern(exper.conditIndexMap,CONST.CONTROL);
    control = AeplUtil.MakeMedianCol(condits(temp));

    exper.t_int = permute(exper.t_int_long(1:size(control.medianCol,1)),[2 1]);

    black = [0,0,0];
    white = [.8,.8,.8];
    red = [1, 0, 0];
    
    laeout = [2,3];
    plotNum = 1;
    figNum = 1;
    
    figure
    
    for group = groups
        if contains(group{1},CONST.CONTROL)
            continue
        end

        gCondits = exper.groupConditMap(group{1});
        
        len = length(gCondits);
        colors_p = [linspace(black(1),white(1),len)', linspace(black(2),white(2),len)', linspace(black(3),white(3),len)'];
    

        subplot(laeout(1),laeout(2),plotNum);
        hold on 

        scatter(exper.t_int,control.medianCol, 10,'filled','MarkerFaceColor',red) ;
        
        subCount = 1;
        for conditNamee = gCondits
            
            conditName = conditNamee{1};
            %disp(conditName)
            if contains(conditName,CONST.CONTROL)
                continue
            end
            
            cond = exper.conditIndexMap(conditName);

            tempCondit = AeplUtil.MakeMedianCol(condits(cond));
            
            scatter(exper.t_int,tempCondit.medianCol, 10,'filled','MarkerFaceColor',colors_p(subCount,:)) ;
            subCount = subCount + 1;
        end
        
        ylim([0 exper.ylimit])
    
        xlabel('time (hours)')             
        ylabel('displacement (pixels)')
        
        temp = [ { CONST.CONTROL }, gCondits ];
        legend(temp)
        
        title(group)
        
        plotNum = plotNum + 1;

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
        
    end
    close(gcf)

end
    
