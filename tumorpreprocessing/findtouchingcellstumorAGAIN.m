load('/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles/merfishMask.mat')%loads mymask
% cellsstats=regionprops(mymask,'Centroid','PixelIdxList','Eccentricity','Area','Circularity','Extent');
% matlabtempfiles='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles'
% numberofrows=size(mymask,1);
% % centroids=zeros(size(cellsstats,1),2);
% % for i=1:size(cellsstats,1)
% %     centroids(i,:)=cellsstats(i).Centroid;
% % end
% % 
% cellsizes=zeros(size(cellsstats,1),2);
% for i=1:size(cellsstats,1)
%     cellsizes(i,:)=cellsstats(i).Area;
% end
% % for i=1:size(cellsizes,1)
% %     centroids(i,:)=[0,0];
% % end
% centroids=table2array(readtable('/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles/coordinates.csv'));
% 
% size(centroids)
% size(cellsstats)
% mydistances=pdist2(centroids,centroids);
% mydistances=0.108*mydistances; 
% mydistances=mydistances<50&mydistances~=0;
% cellsIDs=[1:size(mydistances,2)];
% mydistances=mydistances.*cellsIDs;
% globalFoundSpotsfromBW=cell(size(cellsizes,1),1);
% clear mymask
% 
% save(fullfile(matlabtempfiles,strcat('neighborWorkspace.mat')))
matlabtempfiles='/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles'
numberofrows=size(mymask,1);clear mymask
load(fullfile(matlabtempfiles,strcat('neighborWorkspaceSplit.mat')))
load(fullfile(matlabtempfiles,strcat('neighbordistances.mat')))

delete(gcp('nocreate'))
thatparpool=parpool(19)  
parfor (i=1:(size(cellsizes,1)),19)
% for i=1:(size(cellsizes,1)-1)
    if rem(i,200)==0
        disp(strcat('current cell is ',num2str(i), 'out of ', num2str(size(cellsizes,1))))%progress report
    end
    localneighbors=[];
%     if cellsstats(i).Area>7500
    centercell=cellsstats(i).PixelIdxList;
    centercell=unique([centercell;centercell-30;centercell+30;centercell+30*numberofrows;centercell-30*numberofrows;centercell+30*(numberofrows-1);centercell-30*(numberofrows-1);centercell+30*(numberofrows+1);centercell-30*(numberofrows+1)]);
%     [centercellrow,centercellcol]=ind2sub(size(finalmask),centercell);
    potentialneighbors=mydistances(i,:);
    potentialneighbors=potentialneighbors(potentialneighbors~=0);
    pixelsfromneighbors=[];
    cellitsfrom=[];
      if size(potentialneighbors,2)>0
        for j=1:size(potentialneighbors,2)
            pixelsfromneighbors=[pixelsfromneighbors;cellsstats(potentialneighbors(j)).PixelIdxList];
            cellitsfrom=[cellitsfrom;potentialneighbors(j)*ones(size(cellsstats(potentialneighbors(j)).PixelIdxList,1),1)];
        end
        realneighbors=ismember(pixelsfromneighbors,centercell);
        realneighbors=unique(cellitsfrom(realneighbors));
        localneighbors=[localneighbors;realneighbors];
%     for j=1:size(potentialneighbors,2)
%         checkthatcell=potentialneighbors(j);
%         if cellsstats(checkthatcell).Area>7500
%         thatcell=cellsstats(checkthatcell).PixelIdxList;
%         [thatcellrow,thatcellcol]=ind2sub(size(finalmask),thatcell);
%         distance=pdist2([centercellrow,centercellcol],[thatcellrow,thatcellcol]);
%         if min(min(distance))<50
%             localneighbors=[localneighbors;checkthatcell];
%         end
%         end
%     end
    end
%     end
    globalFoundSpotsfromBW{i}=localneighbors;
end
neighbors=zeros(size(cellsstats,1),size(cellsstats,1));

for i=1:size(cellsstats,1)
    localneighbors=globalFoundSpotsfromBW{i};
    if size(localneighbors,1)>0
        for j=1:size(localneighbors,1)
            neighbors(i,localneighbors(j))=1;
            neighbors(localneighbors(j),i)=1;
        end
    end
end
save('/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles/touchingneighboroverlapUGERAGAIN.mat','neighbors','-v7.3')
numberofneighbors=sum(neighbors,2);
smallertable=zeros(size(neighbors,1),max(max(numberofneighbors)));
for i=1:size(neighbors,1)
    if sum(neighbors(i,:)>0)>0
        nonzeros=find(neighbors(i,:)>0);
        for j=1:sum(neighbors(i,:)>0)
            smallertable(i,j)=nonzeros(j);
        end
    end
end
save('/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles/touchingnumberofneighborsoverlapUGERAGAIN.mat','numberofneighbors','-v7.3')
save('/broad/clearylab/Users/Loic/tumorMERFISHdish1Donor1FFPE/202408041743_tumor1DonorMerfishFFPE_VMSC11302/matlabtempfiles/touchingnumberofneighborsoverlapUGERAGAIN.mat','smallertable','-v7.3')
disp('trying to save csv')
writematrix(smallertable,fullfile(matlabtempfiles,'sample1allcellsTouchingDesignbigcell7500UGERAGAIN.csv'))
writematrix(neighbors,fullfile(matlabtempfiles,'sample1allcellsTouchingDesignbigcell7500UGERAGAIN.csv'))
