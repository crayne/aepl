
deathTest = 'C:\Users\eac84\OneDrive\Desktop\death_ex_for_edgar\death_ex_for_edgar';
LotsOfCells = 'D:\Google\MSK data';

AllDataDrive = 'G:\MSKCC_Baylies-Lab_Rhabdo-Data\*\*\Czi';

[exper,condits] = ReadCsvAsCondits(experPath);

SummarizeData.RunSummarizeData(experPath,exper,condits)
PlotData.RunPlotData(experPath,exper,condits)

ProcessCzi.initConstants
global CONST
CONST.DRAW_MODE = 0;
CONST.DRAW_MODE = 1;


ProcessCzi.RunProcessCzi(LotsOfCells)
profile on 
ProcessCzi.RunProcessCzi(AllDataDrive)
profile viewer

