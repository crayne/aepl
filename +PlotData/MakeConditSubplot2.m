%% MakeConditSubplot  -  beta
%       MakeConditionSubplot
%   given a condition gather data and graph
%
% 


function MakeConditSubplot2(mat,exper,dmso)
    hold on
    
    dmso.mat(dmso.mat>50) = nan;
    dmsoMatRowMeds = nanmedian(dmso.mat,2);
    
    mat(mat>50)=nan;
    %matRowMeans = nanmean(mat,2);
    matRowMeds =  nanmedian(mat,2);
    t_int_temp = permute(exper.t_int(1:size(matRowMeds,1)),[2 1]);
        t_int_temp2 = permute(exper.t_int(1:size(dmsoMatRowMeds,1)),[2 1]);

    %scatter(t_int_temp,matRowMeans, 10,'filled');
    scatter(t_int_temp,matRowMeds, 10,'filled') ;
    scatter(t_int_temp2,dmsoMatRowMeds, 10,'filled') ;

    
    ylim([0 exper.ylimit])
    
    xlabel('time (hours)')              % x-axis label
    ylabel('displacement (pixels)')     % y-axis label
    legend('','dmso')
    
    %l = legend('mean','median');
    %l.Location = 'northeastoutside';
end



