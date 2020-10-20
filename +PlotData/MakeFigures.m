function MakeFigures(Data,well,featureNames)

close all;

treatStr = {'BYL','IKK','PI-103','Tram','Bardox','DMSO .1'};
treatIdx = [2,5,8,11,14,17,21];

%%
Groups = {'RH30','RH41'};
pageNum = 1;



%% Make Tiff Figures
figure('position',[100,100,1800,800]);
filename = 'D:\Google\MSK data\Figures.tif';
delete(filename);
Fr = [];
%% Over Time Figures
if 1

F = PlotData.PlotOverTime(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum);
Fr = [Fr,F];
pageNum = pageNum+length(F);

end
%% Three Panel
if 0 
F = PlotData.PlotOverTime(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum);
Fr = [Fr,F];
pageNum = pageNum+length(F);
end 

if 1 

F = PlotData.PlotStatDiff(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum);
Fr = [Fr,F];
pageNum = pageNum+length(F);
end 

%pgn = ((pageNum-1)/2);
for f = 1:length(Fr)
    
        F = Fr{f};
%         imwrite(F.cdata,filename,'writemode','append');
%         F = Fr{f+pgn};
        imwrite(F.cdata,filename,'writemode','append');
end 



%% Make Figures For the Control 

filename = 'D:\Google\MSK data\FiguresControl.tif';
delete(filename);
Fr = [];

pageNum= 1;
F = PlotData.PlotControlCompare(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum);
Fr = [Fr,F];


for f = 1:length(Fr)
        F = Fr{f};
        imwrite(F.cdata,filename,'writemode','append');
end 

%% Make Death% + Migration Figures For Each Drug 
figure('position',[100,100,1800,800]);
filename = 'D:\Google\MSK data\Figures_DrugBreakDown.tif';
delete(filename);
Fr = [];

pageNum = 1;
F = PlotData.PlotDrugBreakdown(Data,well,featureNames,treatStr,treatIdx,Groups,pageNum);
Fr = [Fr,F];


for f = 1:length(Fr)
        F = Fr{f};
        imwrite(F.cdata,filename,'writemode','append');
end 


