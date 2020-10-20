    

%% constants - for now don't worry about it
global CONST
CONST = struct();

CONST.PLATE_MAP_SUF = '_plate-map.csv';
CONST.CZI_SUF = '.czi';
CONST.CONDIT_DATA_SUF = '-veloc-conditions.csv';
CONST.WELL_DATA_SUF = '-veloc-wells.csv';
CONST.PLOT_SUF = '-veloc-plot.pdf';

CONST.AVG_PREF = 'Avg';
CONST.MEDI_PREF = 'Median';


CONST.CZI_DIR = '';
CONST.CSV_DIR = 'Csv';
CONST.PROCESSED_DIR = 'ProcessedData';
CONST.ANNOTATED_TIF_DIR = 'TifsWithTracks';

%CONST.COL_LAYOUT_VER = '4 cols';
CONST.COL_LAYOUT_VER = '5 cols';
%CONST.COL_LAYOUT = AeplUtil.GetColLayout(CONST.COL_LAYOUT_VER);

CONST.SET_LAEOUT = [2, 3];    % [rows, columns]

%CONST   

CONST.CONTROL = 'DMSO .1%';
CONST.TEMP_PDF = 'temp.pdf';



%% experiment folder names - for now is the folder with all the czi files

% f16_06_23 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/16-06-23 (1st plate)/16-06-23/';
% f17_06_28 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/17-06-28 final pi3k inhibitors/';
% f17_07_13 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/17-07-13 and 20170629 first plate fda screen/7.13.17 Part b/Whitney2/';
% f18_06_20 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/18-06-20/';
% f18_03_18 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/2018-3-18';
% f18_08_02 = '/Volumes/baylieslab/Current Lab Members/Whitney/Whitney  8-2';
% 
% f17_11_10 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/17-11-10 New NF-kB and PI3K FDA approved drugs';
% f17_11_17 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/17-11-17 all cell lines top hits first run molly and whitney';
% f17_12_28 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/17-12-28';


rootFolder = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)';

f18_09_18_name = '20180918/movies';
f18_09_18_folder = fullfile(rootFolder, f18_09_18_name);


f18_09_18 = '/Volumes/baylieslab-1/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH30/20180918/movies';

f18_09_25 = '/Volumes/baylieslab-1/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH30/20180925 (no good)';
f18_09_26 =  '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH30/20180926';
f18_09_27 = '/Volumes/baylieslab-4/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH30/20180927';
f18_09_28 =  '/Volumes/baylieslab-1/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH41/20180928';
f18_09_29 =  '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH41/20180929';
f18_10_03 =  '/Volumes/baylieslab-8/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH41/20181003'
f18_10_23 = '/Volumes/baylieslab-7/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH41/20181023'
f18_10_24 = '/Volumes/baylieslab-7/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH41/20181024/Whitney 10-25'
f18_11_06 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH30/20181106'
f18_11_07 = '/Volumes/baylieslab-1/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH30/20181107'
f18_11_07 = '/Volumes/baylieslab-1/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/Post-mycoplasma data (starting 9:18:18)/RH30/20181109'
%%This is where you put the file name of the czi folder you wan to run. The
%%WTF is an assigned name whereas anything after the = is the path name.
%%For every new one just make a new name ie see below
%fyy_mm_dd = 'folderPathHere';

%% set experPath to the folder you want to run right now from the assigned name above
experPath = f18_11_07;



%% will process the czi and produce folder Cvs with raw data and folder TifsWithTracks with tifs
ProcessCzi.RunProcessCzi(experPath)


%% if folder Csv exists

%% reads in the data
%[exper,condits] = ReadCsvAsCondits(experPath);

%% after data has been read in either create summary csv files or plots
%SummarizeData.RunSummarizeData(experPath,exper,condits)
    %% requires _plate-map.csv
%PlotData.RunPlotData6(experPath,exper,condits)



