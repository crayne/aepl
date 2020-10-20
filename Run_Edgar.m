
deathTest = 'C:\Users\eac84\OneDrive\Desktop\death_ex_for_edgar\death_ex_for_edgar';
LotsOfCells = 'D:\Google\MSK data';

AllDataDrive = 'G:\MSKCC_Baylies-Lab_Rhabdo-Data\*\*\Czi';

TestWell = 'G:\MSKCC_Baylies-Lab_Rhabdo-Data\RH30\18-11-07.rht\Czi';

[exper,condits] = ReadCsvAsCondits(experPath);

SummarizeData.RunSummarizeData(experPath,exper,condits)
PlotData.RunPlotData(experPath,exper,condits)

ProcessCzi.initConstants
global CONST

CONST.DRAW_MODE = 0;
CONST.DRAW_MODE = 1;

ProcessCzi.RunProcessCzi(TestWell)

ProcessCzi.RunProcessCzi(LotsOfCells)

ProcessCzi.RunProcessCzi(AllDataDrive)


RunPlotData();