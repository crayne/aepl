function ExportData(Data,well,featureNames)

filename = 'D:\Google\MSK data\Stats.xlsx';
delete(filename);

 xlswrite(filename,well,'FileList')
for i = 1:size(Data,1)

    for ii = 1
        sheet =[featureNames{i,1} ,' OverTime'];
        A = Data{i,ii};
        %W = 
        xlswrite(filename,A,sheet)
    end
    
end
end