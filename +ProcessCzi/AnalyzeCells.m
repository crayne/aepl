function [cellsOut] = AnalyzeCells(cells)

TIDs = [cells.Tid];
CIDs = [cells.id];

for tid = 1:max(TIDs)
    
  C = cells(TIDs==tid);  
  mitScore = [C.MitoScore];
  mitScore = medfilt1(mitScore,3);
  bMitotic = mitScore>0.07;

  locs = vertcat(C.Centroid);
  dist = [0;sqrt(sum(   (abs(diff(locs)).^2)  ,2))]';
  %% Got Some Signal 
  area = [C.Area]; 
   %TrkPrps = regionprops(bMitotic,'area');
    
cellState = zeros(size(mitScore));  
   
%% Check if Mitotic 
if nnz(bMitotic)>3   

   bMitotic = bMitotic & ~bwareaopen(bMitotic,10);
   bDead = bwareaopen(bMitotic,10); 
  
   cellState(bMitotic) = 1;
   cellState(bDead) = 2;
%%    
end 

%% Check If Dead 
   bDead = dist<5;
   bDead = bwareaopen(bDead,10);
   bDead = bDead | (area < 500);
   bDead = bDead | abs(mitScore)>1;
   
   bDead = bwareaopen(bDead,20);
   
   cellState(bDead) = 2;


%% Check Dead Or Alive Score

if any(cellState)
%% write Mitotic labels 
      ids = [C(cellState == 1).id];
      for i = ids
        cells(CIDs==i).Label = 1;
      end 
      
%% write Mitotic labels  
      ids = [C(cellState == 2).id];
      for i = ids   
        cells(CIDs==i).Label = 2;
      end 
      
end 

end 

cellsOut = cells;
end 


