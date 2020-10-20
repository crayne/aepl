
f18_03_18 = '/Users/baylieslab/Documents/Amelia/rmsMim/18-03-18/';
f18_06_20 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/18-06-20/';
f16_06_23 = '/Volumes/baylieslab/Current Lab Members/Whitney/Rhabdomyosarcoma plate movies/16-06-23 (1st plate)/16-06-23/';


mypath = f16_06_23;
%mypath = f18_06_20;


tic


% [exper,condits] = readd(mypath);
%[exper,condits] = readd(f18_06_20);

exper.frames = length(condits(1).wells(1).cells(1).xcoords);
[exper,condits] = plotCondits(exper, condits);


toc
