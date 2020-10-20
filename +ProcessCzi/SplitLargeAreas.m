%% 
function imBWout = SplitLargeAreas(BW)

%     AreaT = round(numel(BW)*0.004);
%     AreaT2 = round(numel(BW)*0.002);

%% Size threshhold to apply split 
    AreaT3 = 1500;
    AreaT2 = 400; %% First Threshold for Watershed < BEST option 
    AreaT1 = 600;

    imBW = imresize(BW,0.5);
    
    %% Split Step 1 
    imBWout = uint8(imBW);
    tooLarge = bwareaopen(imBW, AreaT2);
    imBWout(tooLarge) = 0;
    
    bI = ProcessCzi.SeparateObjects1(tooLarge);
    
    %% Split Step 2 
    tooLarge = bwareaopen(bI, AreaT3);
    bPassed = bI & ~tooLarge;
    imBWout(bPassed) = 2;
    
    L = ProcessCzi.SeparateObjects(tooLarge,AreaT1);

    imBWout(L>0) = 2;
   
%     figure 
%     imagesc([imBW,imBWout, bI])
%     axis equal 
%     figure 
%     imagesc([label2rgb(bwlabel(imBW),'jet','k','shuffle'),label2rgb(bwlabel(bI),'jet','k','shuffle')])
%     
    imBWout = imresize(imBWout,size(BW),'method','nearest');
    
    
end
