%% DrawTracks  -  beta
%   from Epl
%   I've made a fair amount of changes
%

function DrawTracks(im,Segs,fName,videoHandle)
fprintf(1,'\tDrawing tracks and saving to tif\n')
fprintf(1,'\t\tframe 1')



cmap = colormap('jet');
cmap = cmap(randperm(size(cmap,1)),:);

saveframeTiff  = 0;
if exist('videoHandle','var') && ~isempty(videoHandle)
    saveframeMP4  = 1;
else
    saveframeMP4  = 0;
end

AllTracks = [Segs.Tid];
Tracks = unique(AllTracks);
counts = hist(AllTracks,Tracks);

%%
[~,fnameShort,~]  = fileparts(fName);


Fr = {};
h = figure('position',[50 50 1400 700]);

if 0
    % plot
    
    imt = im(:,:,1);
    imt = imadjust(mat2gray(imt));
    
    imagesc([imt,imt])
    colormap gray
    hold on
    T = [Segs.Tid];
    for i = 1:max(T)
        
        if counts(Tracks==i)<5
            continue
        end
        bT = T==i;
        if nnz(bT)==0; continue; end
        
        cid = mod(i,64)+1;
        
        trackcolor = cmap(cid,:);
        pts = vertcat(Segs(bT).Centroid);
        plot(pts(:,1),pts(:,2),'-','linewidth',2,'color',trackcolor )
        
    end
    
    F = getframe(h);
    F = F.cdata;
    F = imresize(F,0.5);
    Fr{end+1} = F;
    
    clf()
    
    ax1 = axes('Position',[0.05 0.05 0.45 0.9]);
    ylabel('Time')
    title('X')
    hold on
    ax2 = axes('Position',[0.55 0.05 0.45 0.9]);
    title('Y')
    hold on
    T = [Segs.Tid];
    for i = 1:max(T)
        
        if counts(Tracks==i)<1
            continue
        end
        
        cid = mod(i,64)+1;
        
        trackcolor = cmap(cid,:);
        
        
        bT = T==i;
        if nnz(bT)==0; continue; end
        pts = vertcat(Segs(bT).Centroid);
        times = vertcat(Segs(bT).time);
        
        
        plot(ax1,pts(:,1),times,'-','linewidth',2,'color',trackcolor )
        plot(ax2,-pts(:,2),times,'-','linewidth',2,'color',trackcolor )
        
    end
    
    F = getframe(h);
    F = F.cdata;
    Fr{end+1} = F;
    
end





%%
%h = figure('Visible', 'off');
% h = figure('position',[0 0 1600 700]);
for i = 1:size(im,3)
    clf(h,'reset')
    %set(h, 'Visible', 'off')
    %h = figure('Visible', 'off');
    %h = figure();
    
    Tsegs = Segs([Segs.time]==i);
    % put break on next line to step through outlined cells
    TTracks = [Tsegs.Tid];
    
    imt = im(:,:,i);
    imt = imadjust(mat2gray(imt));
    
    imagesc([imt,imt])
    colormap gray
    Axis = gca;
    Axis.Position = [0 ,0, 1, 1];
    axis equal
    axis off
    hold on
    
    for ii = 1:length(TTracks)
        
        pts = vertcat(Tsegs(ii).Bound);
        cid = mod(TTracks(ii),64)+1;
        
        cent = Tsegs(ii).Centroid;
        
        xs = size(imt,1);
        cent(1) = cent(1)+xs;
        
        if counts(Tracks==Tsegs(ii).Tid)>10
            trackcolor = cmap(cid,:);
        else
            trackcolor = [0.5,0.5,0.5];
        end
        
        plot(pts(:,2)+xs,pts(:,1),'-','Color',trackcolor,'LineWidth',2)
        textcolor = [1,1,1]
        text(cent(1),cent(2)-5,num2str(Tsegs(ii).Tid),'Color',textcolor,'FontSize',10,'BackgroundColor', [.5,.5,.5])
        
        if Tsegs(ii).Label==1
            plot(cent(1),cent(2),'+','Color','g','markersize', 20)
%             viscircles([cent(1),cent(2)],20,'Color','g');
        elseif Tsegs(ii).Label==2
            plot(cent(1),cent(2),'x','Color','r','markersize', 20)
%             viscircles([cent(1),cent(2)],20,'Color','r');
        end
    end
    
    text(10,10,fnameShort,'color','red')
    %hold off
    
    
    
    drawnow
    F = getframe(h);
    F = F.cdata;
    Fr{end+1} = F;
    
end
close(h)

%%
if saveframeTiff
    firstRound = 1;
    for i = 1:length(Fr)
        if firstRound == 1
            imwrite(Fr{i},fName,'WriteMode','overwrite')
            firstRound = 0;
        else
            imwrite(Fr{i},fName,'WriteMode','append')
        end
        
        %
    end
end
%%

if saveframeMP4
    
    for i = 1:length(Fr)
        writeVideo(videoHandle,Fr{i})    
    end
end

end
