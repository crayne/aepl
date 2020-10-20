function [exper,condits2] = plotCondits(exper,condits)
    SUB_DIR = 'pdfPlots';
    PLOT_SUF = '-plot.pdf';

    myPath = fullfile(exper.folder, SUB_DIR);
    if ~exist(myPath,'dir')
        mkdir(myPath)
    end    

    exper.t_int = (1/6):(1/6):(exper.frames+10/6);

    exper.ylimit = 20;
    
    groups = exper.groupWellMap.keys();
   
    condits2(length(condits)) = struct('name','','wells',[],'mat',[]);

    for group = groups
        gCondits = exper.groupWellMap(group{1});
        plotCount = length(gCondits);
        laeout = makeLaeout(plotCount);
        
        figure
        plotNum = 1;
        for conditNamee = gCondits
            conditName = conditNamee{1};
            cond = exper.conditIndexMap(conditName);
%              try
               temp  = makeConditMat(condits(cond));
               condits2(cond)= temp;
%             catch e
%                 disp(conditName)
%                 disp(e.getReport)
%             end
            subplot(laeout(1),laeout(2), plotNum)
            title(conditName)
            
            plottyPlot(condits2(cond).mat,exper);
            plotNum = plotNum + 1;
        end
%         h = gcf;
%         set(h, 'PaperUnits','inches','PaperPosition',[0 0 11 8.5],'PaperOrientation','landscape')
%         
%         
%         fname = [exper.name,'_',group{1},PLOT_SUF];
%         
%         print(h,fullfile(myPath,fname),'-dpdf')
%         close(h)
    end
    %condits = condits2;
end


function laeout = makeLaeout(plotCount)
    if plotCount == 1
       laeout = [3,4]; 
   
    else
        colCount = floor(sqrt(plotCount));
        rowCount = colCount;
        n = plotCount - colCount*rowCount;
        if n > 0 
            colCount = colCount + ceil(n/rowCount);
        end
        laeout = [rowCount colCount];
    end
end


function condit = makeConditMat(condit)
    condit.mat = [];
    col = 1;
    for w = 1:length(condit.wells)
        try 
            frames = length(condit.wells(w).cells(1).distance(2:end));

            condit.wells(w).mat = [];
            wCol = 1;

            for c = 1:length(condit.wells(w).cells)
                try 
                    condit.mat(1:frames,col) = condit.wells(w).cells(c).distance(2:end);
                    col = col+1;


                    condit.wells(w).mat(1:frames,wCol) = condit.wells(w).cells(c).distance(2:end);
                    wCol = wCol + 1;

                    condit.wells(w).cells(c).coords = horzcat(condit.wells(w).cells(c).xcoords, condit.wells(w).cells(c).ycoords);


                catch e
                    fprintf(2,['condition: ', condit.name,'\n'])
                    fprintf(2,[e.getReport(),'\n'])
                end
            end
        catch e
            fprintf(2,['condition: ', condit.name,'\n'])
            fprintf(2,[e.getReport(),'\n'])
        end
        
        
        
        
    end
end


function plottyPlot(mat,exper)
    hold on
    
    mat(mat>50)=nan;
    matRowMeans = nanmean(mat,2);
    matRowMeds =  nanmedian(mat,2);
    t_int_temp = permute(exper.t_int(1:size(matRowMeans,1)),[2 1]);
    
    scatter(t_int_temp,matRowMeans, 10,'filled');
    scatter(t_int_temp,matRowMeds, 10,'filled') ;
    
    ylim([0 exper.ylimit])
    
    xlabel('time (hours)') % x-axis label
    ylabel('displacement (pixels)') % y-axis label
    legend('mean','median')
    %l = legend('mean','median');
    %l.Location = 'northeastoutside';
end
