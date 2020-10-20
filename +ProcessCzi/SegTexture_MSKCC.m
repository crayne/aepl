%%  SegTexture_MSKCC  -  Epl
%       from Epl
%       little to no changes on my part
%


function [CC,imBW] = SegTexture_MSKCC(Iin,t,DrawPlot)
if isempty(Iin)
    CC = [];
    imBW = [];
    return
end

% Disk7 = strel('disk',7);
% Disk3 = strel('disk',3);
Disk20 = true(15);
Disk11 = true(11);
Disk5 = true(5);

MaxArea = 1800;
MinArea = 500;
MinArea2 = 50;

%% Filter Image
I = single(Iin);
I = I - median(I(1:2:end));
% I(~imdilate(imBW,Disk7))=0;
I = I/std(I(:));
I = min(I,3);
I = max(I,-3);
I = I/std(I(:));

%% Filter using FFT ??
% fI = (fft2(I));
% mag = abs(fft2(I));
% mag1 = fftshift(mag);
% mag1 = imerode(mag1,ones(3));
% fI = fI.* mag1;
% I2 = ifft2(fI);


%% Transform using a texture filter
I1 = mat2gray(I);

Igrad1 = I1;

Igrad1 = imresize(Igrad1,0.5);
Igrad1 = single(Igrad1);
Igrad1 = stdfilt(Igrad1);

Igrad1 = single(Igrad1);
Igrad1 = imresize(Igrad1,size(I1));
Igrad1 = imadjust(mat2gray(Igrad1));

%I2 = single(Igrad1);
I2 = mat2gray(Igrad1 + mat2gray(I));
%% Segment Cells
[h,counts] = histcounts(I2,100);

% cumh = cumsum(h);
% cumh = cumh / cumh(end);
% testHigh = counts(find(cumh>0.85,1));
i = 1:length(h);
maxH = max(h);
maxidx = find(h == maxH,1);
tidx = find( h < maxH*0.1 &  i >maxidx , 1);

testHigh = counts( tidx );
thresh = [0,testHigh];

%% Morphological operations
imBW = I2 >= thresh(2);
imBW = bwareaopen(imBW,5);
imBW = imclose(imBW,Disk5);
imBW1 = imopen(imBW,Disk5);
imBW = imerode(imBW1,Disk20);
imBW = bwareaopen(imBW,MinArea);


ja=1;
subaxis('Spacing',0,'Padding',0, 'Margin',0)

subaxis(2,3,ja)
imagesc(I)
ja=ja+1;

subaxis(2,3,ja)
imagesc(imBW)
ja=ja+1;

%% Find the Edges that separate the cells
bEdges = imdilate(imBW,ones(5)) & ~imerode(imBW,ones(5));
EdgeVal = median(I(bEdges));

subaxis(2,3,ja)
imagesc(bEdges)
ja=ja+1;




if EdgeVal>1
    bD = I>1;
else
    bD = I<-1;
end

bD = bwareaopen(bD,20);

subaxis(2,3,ja)
imagesc(bD)
ja=ja+1;


% bD = bD & ~imopen(imclose(bD,Disk7),Disk3);
% bL = I>4;
% bL = bwareaopen(bL,MinArea2);
% bL = bL & ~imopen(bL,Disk10);
bE = bD;
bE = imopen(bE,ones(1));
bE = bwareaopen(bE,MinArea2);
%bE = imdilate(bE,strel('disk',1));

%% Mitotic Segmentation
bBright = I>2;
bBright(~imBW1) = 0;

bBright = bwareaopen(bBright,100);
bBright = imclose(bBright,Disk5);
bBright = imopen(bBright,Disk11);

bBright = bwareaopen(bBright,100);
bBright = imdilate(bBright,Disk5);

bBrtE = imdilate(bBright,Disk11)& ~bBright;

%% Remove to make final Binary image

bE = bE & ~bBright;

imBW = imBW & ~bBrtE;

imBWbig = bwareaopen(imBW,1800);
imBWsml = imBW & ~imBWbig;
imBW = imBWsml | (imBWbig & ~bE);
imBW = imfill(imBW,'holes');

imBW = bwareaopen(imBW,MinArea);

imBWFinal = imBW | bBright;
imBWFinal(1:50,1:50) = 0;

%% Split Segmentation
L1 = ProcessCzi.SplitLargeAreas(imBWFinal);
imBW = L1>0;
imBW = bwareaopen(imBW,200);

CC = ProcessCzi.MakeCells(imBW,t, double(bBright), double(L1));





%%  Show Segmentation
if DrawPlot
    
    I3 = cat(3,imBW,bE,bBright);
    %ProcessCzi.Plot4(I2,testHigh,I,I1,I3)
    
    
    ProcessCzi.Plot6(Igrad1,testHigh,I,I1,I3)
    ax = subplot(2,3,6);
    L1(~imBW) = 0;
    LC = mat2gray(label2rgb(bwlabeln(L1),'jet','k','shuffle'));
    Icolor = ProcessCzi.MakeColor(Iin,LC);
    
    imagesc(ax,Icolor)
    ax.Position = [2/3,0,1/3,0.5];
    grid on
    
end

end

