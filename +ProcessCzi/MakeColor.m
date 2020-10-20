

function Icolor = MakeColor(I1,bI3)

bI3 = im2uint8(bI3);

bR = bI3(:,:,1);
bG = bI3(:,:,2);
bB = bI3(:,:,3);

Icolor = im2uint8(imadjust(mat2gray(I1)));
Icolor = repmat(Icolor,[1 1 3]);
Icolor = Icolor*0.7 + 0.2/255;


color = [1 -1 -1];
for i = 1:3
    Icolor(:,:,i) = Icolor(:,:,i) +  0.5*bR*color(i);
end

color = [-1 1 -1];
for i = 1:3
    Icolor(:,:,i) = Icolor(:,:,i) +  0.8*bG*color(i);
end

color = [-1 -1 1];
for i = 1:3
    Icolor(:,:,i) = Icolor(:,:,i) +  0.5*bB*color(i);
end
end