load('/broad/clearylab/Users/Loic/thp1homemadezombie_1/matlabtempfiles/finalmask.mat')%load cell mask
cellsstats=regionprops(finalmask,'Centroid','PixelIdxList','Eccentricity','Area','Circularity','Extent');%load list of cells and their stats
matlabtempfiles='/broad/clearylab/Users/Loic/thp1homemadezombie_1/matlabtempfiles'

centroids=zeros(size(cellsstats,1),2);%get centroid to first find cells that have close centroids
for i=1:size(cellsstats,1)
    centroids(i,:)=cellsstats(i).Centroid;
end

cellsizes=zeros(size(cellsstats,1),2);
for i=1:size(cellsstats,1)l
    cellsizes(i,:)=cellsstats(i).Area;
end
for i=1:size(cellsizes,1)%remove cells that are too small (like we did when filtering the cout table , these are debris
    if cellsizes(i)<7500
        centroids(i,:)=[0,0];
    end
end
distances=pdist2(centroids,centroids);
distances=0.108*distances; %convert pixels to um
distances=distances<50&distances~=0;%keep cells that are closeby, remove identity
cellsIDs=[1:size(distances,2)];
distances=distances.*cellsIDs;
globalFoundSpotsfromBW=cell(size(cellsizes,1),1);
delete(gcp('nocreate'))
thatparpool=parpool(50)  
parfor (i=1:(size(cellsizes,1)-1),22)%in a parallel loop, go over each putative neighbor, and see if their masks actually touch, this did significantly improve the signal
    if rem(i,100)==0
        disp(strcat('current Cell is ',num2str(i)))%progress report
    end
    localneighbors=[];
    if cellsstats(i).Area>7500
        centercell=cellsstats(i).PixelIdxList;
        [centercellrow,centercellcol]=ind2sub(size(finalmask),centercell);%convert index positin back to x,y coordinates for distance computation
        potentialneighbors=distances(i,:);
        potentialneighbors=potentialneighbors(potentialneighbors~=0);
        if size(potentialneighbors,2)>0%if the cell has at least on putative neighbor
            for j=1:size(potentialneighbors,2)%we check them one at a time
                checkthatcell=potentialneighbors(j);
                if cellsstats(checkthatcell).Area>7500
                    thatcell=cellsstats(checkthatcell).PixelIdxList;
                    [thatcellrow,thatcellcol]=ind2sub(size(finalmask),thatcell);
                    distance=pdist2([centercellrow,centercellcol],[thatcellrow,thatcellcol]);%compute distance from each pixels of one cell to each pixel of the other
                    if min(min(distance))<50
                        localneighbors=[localneighbors;checkthatcell];%if cells have at least one pair pof pixels less that 50 pixels away, they touch
                    end
                end
            end
        end
    end
    globalFoundSpotsfromBW{i}=localneighbors;
end
neighbors=zeros(size(cellsstats,1),size(cellsstats,1));

for i=1:size(cellstats,1)
    localneighbors=globalFoundSpotsfromBW{i};
    if size(localneighbors,1)>0
        for j=1:size(localneighbors,1)
            neighbors(i,localneighbors(j))=1;
            neighbors(localneighbors(j),i)=1;
        end
    end
end
save('/broad/clearylab/Users/Loic/thp1homemadezombie_1/properneighbors/touchingneighborcompleteparallel.mat','neighbors','-v7.3')
numberofneighbors=sum(neighbors,2);
save('/broad/clearylab/Users/Loic/thp1homemadezombie_1/properneighbors/touchingnumberofneighborsparallel.mat','numberofneighbors','-v7.3')
disp('trying to save csv')
writematrix(neighbors,fullfile(matlabtempfiles,'sample1allcellsTouchingDesignbigcell7500.csv'))

% %% increase distance to 80 pixels
% disp('increase distance')
% neighbors=zeros(size(cellsstats,1),size(cellsstats,1));
% 
% centroids=zeros(size(cellsstats,1),2);
% for i=1:size(cellsstats,1)
%     centroids(i,:)=cellsstats(i).Centroid;
% end
% 
% cellsizes=zeros(size(cellsstats,1),2);
% for i=1:size(cellsstats,1)
%     cellsizes(i,:)=cellsstats(i).Area;
% end
% for i=1:size(cellsizes,1)
%     if cellsizes(i)<3000
%         centroids(i,:)=[0,0];
%     end
% end
% distances=pdist2(centroids,centroids);
% distances=0.108*distances; 
% distances=distances<50&distances~=0;
% cellsIDs=[1:size(distances,2)];
% distances=distances.*cellsIDs;
% 
% for i=1:(size(cellsstats,1)-1)
%     if rem(i,100)==0
%         disp(strcat('current Cell is ',num2str(i)))
%     end
%     if cellsstats(i).Area>7500
%     centercell=cellsstats(i).PixelIdxList;
%     [centercellrow,centercellcol]=ind2sub(size(finalmask),centercell);
%     potentialneighbors=distances(i,:);
%     potentialneighbors=potentialneighbors(potentialneighbors~=0);
%     if size(potentialneighbors,2)>0
%     for j=1:size(potentialneighbors,2)
%         checkthatcell=potentialneighbors(j);
%         if cellsstats(checkthatcell).Area>7500
%         thatcell=cellsstats(checkthatcell).PixelIdxList;
%         [thatcellrow,thatcellcol]=ind2sub(size(finalmask),thatcell);
%         distance=pdist2([centercellrow,centercellcol],[thatcellrow,thatcellcol]);
%         if min(min(distance))<50
%             neighbors(i,checkthatcell)=1;
%             neighbors(j,checkthatcell)=1;
%         end
%         end
%     end
%     end
%     end
% end
% writematrix(neighbors,fullfile(matlabtempfiles,'sample1allcellsTouchingDesignbigcell7500_80pixels.csv'))
% save('/broad/clearylab/Users/Loic/thp1homemadezombie_1/properneighbors/touchingneighborcomplete80pixels.mat','neighbors','-v7.3')
% numberofneighbors=sum(neighbors,2);
% save('/broad/clearylab/Users/Loic/thp1homemadezombie_1/properneighbors/touchingnumberofneighbors80pixels.mat','numberofneighbors','-v7.3')
