x = -10:.1:10;
y1 = sin(x);
y2 = cos(x);
clf
hold on

%f1 = figure('Visible','off');
%ax1 = axes(f1,'Visible','off');


p1 = plot(x,y1);
p2 = plot(x,y2);
% p1 = plot(ax1,x,y,'Visible','off');
