
function bOut = SeparateObjects1(imBW)

if nnz(imBW)==0
    bOut = imBW;
    return
end

imBW2 = imerode(imBW,strel('disk',3));
dist = bwdist(~imBW2);
BW2 = imregionalmax(dist);
[~,I] = bwdist(BW2);

maxD = dist;
maxD(imBW) = dist(I(imBW));

Erode = dist > maxD*0.5;
Erode = bwareaopen(Erode,10);

D = bwdistgeodesic(imBW, Erode);
D(~imBW) = Inf;
L3 = watershed(D,8);
L3(~imBW) = 0;

bOut = L3>0;

%imagesc([imBW,Erode,bOut])


end