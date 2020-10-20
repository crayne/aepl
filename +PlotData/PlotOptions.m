%% PlotOptions  -  alpha
%   a partially constant structure
%   information about different plots
%       descrip (discription) -> constant
%       make -> boolean of whether or not to make that plot type
%



function plotOptions = PlotOptions
    plotOptions = struct();
    
    plotOptions.laeoutByGroup.make = false;
    plotOptions.laeoutByGroup.descrip = [...
        'Rows and columns for each group are calculated given the number of plots in that group\n',...
        'Creates a shape close to a square'...
    ];
    
    plotOptions.laeoutByOneGroup.make = false;
    plotOptions.laeoutByOneGroup.descrip = [...
        'Rows and columns are caluclated given the number of plots in the largest group\n',...
        'Creates a shape close to a square'...
    ];

    plotOptions.laeoutSet.make = true;
    plotOptions.laeoutSet.descrip = [...
        'Uses CONST.SET_LAYOUT for all groups\n',...
        'Creates multiple figures for groups with more plots than the laeout total'...
    ];
    
end

