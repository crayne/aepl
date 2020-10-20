%% ExportTrackStats  -  beta
%   from Epl
%   I've made a fair amount of changes
%


function ExportTrackStatsSmall(Cells,Dims,filePath)
fprintf(1,'\tMaking csv files\n')

writesql = 0;

%%
%     [P,F,~] = fileparts(Name);
%     P2 = fullfile(P,'Xls');
%     if ~exist(P2,'dir')
%         mkdir(P2)
%     end
%     XlsName = fullfile(P2,[F,'.csv']);

Mat = {};
numcell = size(Cells,1);

Mat(1,1) = {'Cell ID'};
Mat(1,2) = {'x(Pixel Position)'};
Mat(1,3) = {'y(Pixel Position)'};
Mat(1,4) = {'Area(pixels^2)'};
Mat(1,5) = {'Time(Frame Num)'};
Mat(1,6) = {'Track ID'};
Mat(1,7) = {'Label'};
Mat(1,8) = {'wasSplit'};
Mat(1,9) = {'EdgeCost'};

Cents = vertcat(Cells.Centroid);

Mat(2:numcell+1,1) = num2cell([Cells.id]);
Mat(2:numcell+1,2) = num2cell(Cents(:,1));
Mat(2:numcell+1,3) = num2cell(Cents(:,2));
Mat(2:numcell+1,4) = num2cell([Cells.Area]);
Mat(2:numcell+1,5) = num2cell([Cells.time]);
Mat(2:numcell+1,6) = num2cell([Cells.Tid]);
Mat(2:numcell+1,7) = num2cell([Cells.Label]);
Mat(2:numcell+1,8) = num2cell([Cells.wasSplit]);
Mat(2:numcell+1,9) = num2cell([Cells.Edge]);
%%
T = cell2table(Mat);

writetable(T,filePath,'WriteVariableNames',0);

%%

if writesql
    
    [p,f,~] = fileparts(filePath);
    sqlpath = fullfile(p,[f,'.sqlite']);
    
    conn = database(sqlpath, '','', 'org.sqlite.JDBC', 'jdbc:sqlite:');
    
    VNames = regexprep(Mat(1,:),'[ ()^]','_');
    T = cell2table(Mat(2:end,:),'VariableNames',VNames);
    sqlwrite(conn,'Cells',T)
    
    close(conn)
end
%disp(Mat)
%xlswrite(XlsName,Mat)
end