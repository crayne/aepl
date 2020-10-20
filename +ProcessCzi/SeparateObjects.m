
function L2 = SeparateObjects(imBW,MeanArea)


L = bwlabel(imBW);
L2 = zeros(size(L));
if nnz(imBW)==0
    return 
end 

CellFeats = regionprops(imBW,'Area','PixelList','PixelIdxList');
LMat = zeros(size(imBW,1),size(imBW,2),size(imBW,3));
Lmax = 1;


N = 3;
Erode = imerode(imBW,ones([N,N,1]));

Cent = [];
for i = 1:max(L(:))
    k =  round((CellFeats(i).Area/MeanArea));
    if k < 2
        L2(CellFeats(i).PixelIdxList)=Lmax ;
        LMat = LMat + bwperim(L2==Lmax);
        Lmax = Lmax + 1;
    else
        
        idxlist = CellFeats(i).PixelIdxList;

        
        temp = false(size(imBW,1),size(imBW,2),size(imBW,3));
        temp(idxlist) = 1;
        T1 = temp;
        T1 = T1&Erode;

        N = 10;
        [r,c] = find(temp);
        xyList = [r,c];
        xyList = xyList(1:N:end,:);
        
        %%
        try
            T = kmeans(xyList,k);
        catch e
            fprintf(2,'kmeans error, I dont know whats causing it')
            fprintf(2,"exception: " + getReport(e)+"\n")
        end
        
        Tcent = zeros(max(T(:)),2);
        %%
        for j = 1:max(T(:))
        Tcent(j,1) = mean(xyList(T==j,1));
        Tcent(j,2) = mean(xyList(T==j,2));
        end 
        
        Cent = [Cent;Tcent];
    end  
end 
        
        D = bwdistgeodesic(imBW, round(Cent(:,2)), round(Cent(:,1)));
        D(~imBW) = Inf;
        L3 = watershed(D,8);
        L3(~imBW) = 0;
        
        %[~,I] = pdist2(Tcent(:,[2,1]),CellFeats(i).PixelList,'euclidean','Smallest',1);
        
% %         for j = 1:max(L3(:))
% %             
% %             L2(L3==j)=Lmax ;
% %             LMat = LMat + imdilate(bwperim(L2==Lmax),ones(3,3,1));
% %             Lmax = Lmax + 1;
% %         end
        
        %%
%         imagesc(L2);
%         hold on 
%         plot(Tcent(:,2),Tcent(:,1),'ko')

L2 = L3;

% %% Seperate Touching Objects
% Boundry = imdilate(LMat>1,ones([3,3,1]));
% L2(Boundry) = false;
end 