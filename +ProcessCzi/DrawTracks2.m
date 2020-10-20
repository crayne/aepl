%% DrawTracks  -  beta
%   from Epl
%   I've made a fair amount of changes
%

function DrawTracks2(im,Segs,fName)
fprintf(1,'\tDrawing tracks and saving to tif 2\n')
ugh = 0;


close all, 
figure, 

cmap = colormap('jet');
cmap = cmap(randperm(size(cmap,1)),:);

saveframe  = 1;
AllSegs = vertcat(Segs{:});

AllTracks = [AllSegs.Tid];
Tracks = unique(AllTracks);
% counts = hist(AllTracks,Tracks);

Fr = cell(1,size(im,3));


imagesc(im(:,:,1))
    colormap gray
    Axis = gca;
    Axis.Position = [0 ,0, 1, 1];
    axis equal
    axis off
    hold on
    
    
for i = 1:10:size(im,3)
    
    Tsegs = AllSegs([AllSegs.time]==i);
    % put break on next line to step through outlined cells
    TTracks = [Tsegs.Tid];
    
    
    for ii = 1:length(TTracks)
        
        pts = vertcat(Tsegs(ii).Bound);
        cid = mod(TTracks(ii),64)+1;
        plot(pts(:,2),pts(:,1),'-','Color',cmap(cid,:))
        cent = Tsegs(ii).Centroid;
        text(cent(1),cent(2),num2str(Tsegs(ii).Tid),'Color',cmap(cid,:),'FontSize',10)
        
    end
    %hold off
    drawnow
    
    
    %fName = 
  
end
  if saveframe
        %Fr{i} = getframe(gca);
        tempGca = getframe(gca);
        if ugh == 0
            imwrite(tempGca.cdata,fName)
            ugh = 1;
        else
            imwrite(tempGca.cdata,fName,'WriteMode','append') 
        end
    end
%MakeMovie(Fr,Name)
end
