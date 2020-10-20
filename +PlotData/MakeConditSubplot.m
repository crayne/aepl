%% MakeConditSubplot  -  beta
%       MakeConditionSubplot
%   given a condition gather data and graph
%
% 


function MakeConditSubplot(mat,exper)
    hold on
    
    mat(mat>50)=nan;
    matRowMeans = nanmean(mat,2);
    matRowMeds =  nanmedian(mat,2);
    t_int_temp = permute(exper.t_int(1:size(matRowMeans,1)),[2 1]);
    
    scatter(t_int_temp,matRowMeans, 10,'filled');
    scatter(t_int_temp,matRowMeds, 10,'filled') ;
    
    ylim([0 exper.ylimit])
    
    xlabel('time (hours)')              % x-axis label
    ylabel('displacement (pixels)')     % y-axis label
    legend('mean','median')
    
    %l = legend('mean','median');
    %l.Location = 'northeastoutside';
end



