function Plot6(Igrad1,testHigh,I,I1,I3)

I = imresize(I,0.5);
I1 = imresize(I1,0.5);
I3 = imresize(I3,0.5);
Igrad1 = imresize(Igrad1,0.5);

subplot(4,3,2)
histogram(Igrad1(1:10:end))
hold on
plot([testHigh,testHigh], ylim)
hold off


subplot(4,3,5)
histogram(I(1:5:end))
hold on
plot([1,1], ylim)
plot([-1,-1], ylim)
xlim([-3 3])
grid on
hold off
%
subplot(2,3,3)
ID = -I;
ID(ID<0) = 0;
IL = I;
IL(IL<0) = 0;
Icolor = cat(3,mat2gray(IL),mat2gray(ID),I*0);
imagesc(Icolor)
ax = gca();
ax.Position = [2/3,0.5,1/3,0.5];

subplot(2,3,1)
imagesc(I1)
colormap gray
grid on
ax = gca();
ax.Position = [0,0.5,1/3,0.5];

subplot(2,3,4)
imagesc(Igrad1)
colormap gray
grid on
ax = gca();
ax.Position = [0,0,1/3,0.5];

subplot(2,3,5)
Icolor = ProcessCzi.MakeColor(I1,I3);
imagesc(Icolor)
grid on 
ax = gca();
ax.Position = [1/3,0,1/3,0.5];
