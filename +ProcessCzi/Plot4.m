function Plot4(Igrad1,testHigh,I,I1,I3)

I = imresize(I,0.5);
I1 = imresize(I1,0.5);
I3 = imresize(I3,0.5);
Igrad1 = imresize(Igrad1,0.5);

subplot(2,2,2)
ID = -I;
ID(ID<0) = 0;
IL = I;
IL(IL<0) = 0;
Icolor = cat(3,mat2gray(IL),mat2gray(ID),I*0);
imagesc(Icolor)
ax = gca();
ax.Position = [1/2,0.5,1/2,0.5];

subplot(2,2,1)
imagesc(I1)
colormap gray
grid on
ax = gca();
ax.Position = [0,0.5,1/2,0.5];

subplot(2,2,3)
imagesc(Igrad1)
colormap gray
grid on
ax = gca();
ax.Position = [0,0,1/2,0.5];

subplot(2,2,4)
Icolor = ProcessCzi.MakeColor(I1,I3);
imagesc(Icolor)
grid on 
ax = gca();
ax.Position = [1/2,0,1/2,0.5];
