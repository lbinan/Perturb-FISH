mypath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/data'
matlabtempfiles='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/matlabtempfiles'
nucleipath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/nuclei'
load(fullfile(matlabtempfiles,'mymaskregion0V2.mat'));
load(fullfile(matlabtempfiles,'globalmerfishspotsregion0.mat'));
globalmerfishspots=globalmerfishspots(globalmerfishspots(:,end-3)<112000,:);
globalmerfishspots=globalmerfishspots(globalmerfishspots(:,end-2)<71190,:);
% for i=1:size(globalmerfishspots,1)
%     spotsimage(floor(globalmerfishspots(i,end-3)),floor(globalmerfishspots(i,end-2)))=globalmerfishspots(i,1);
% end
mymask=mymask>0;
%  ids=find(spotsimage)>0;
ids=sub2ind([size(mymask,1),size(mymask,2)],globalmerfishspots(:,end-3)-41600+1,globalmerfishspots(:,end-2)-12206+1);
globalmerfishspots(:,13)=ids;
cellstats=regionprops(mymask,'Area','PixelIdxList');
for i=1:size(cellstats,1)
    mypixels=cellstats(i).PixelIdxList;
    positions=ismember(globalmerfishspots(:,13),mypixels);
    globalmerfishspots(positions,14)=i;%add a cell ID to each transcript
end
countable=zeros(size(cellstats,1),3+max(globalmerfishspots(:,2)));
for thiscell=1:size(cellstats,1)%build the count table
    if cellstats(thiscell).Area>0 
        countable(thiscell,1)=thiscell;%#index of that cell
        countable(thiscell,2)=cellstats(thiscell).Area;
        mytable2=globalmerfishspots(globalmerfishspots(:,14)==thiscell,:);
        for gene=1:max(globalmerfishspots(:,2))+1
            cellcountable(1,gene)=sum(mytable2(:,2)==gene-1);
        end
        countable(thiscell,3:end)=cellcountable;
    end
end
writematrix(countable,fullfile(matlabtempfiles,'merfishcounttableregion0rightIndex.csv'))
%%
clear
mypath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/data'
matlabtempfiles='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/matlabtempfiles'
nucleipath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/nuclei'
load(fullfile(matlabtempfiles,'mymaskregion1V2.mat'));
load(fullfile(matlabtempfiles,'globalmerfishspotsregion1.mat'));
globalmerfishspots=globalmerfishspots(globalmerfishspots(:,end-3)<116000,:);
globalmerfishspots=globalmerfishspots(globalmerfishspots(:,end-2)<54400,:);
% spotsimage=zeros(size(mymask));
% for i=1:size(globalmerfishspots,1)
%     spotsimage(floor(globalmerfishspots(i,end-3)),floor(globalmerfishspots(i,end-2)))=globalmerfishspots(i,1);
% end
mymask=mymask>0;
%  ids=find(spotsimage)>0;
ids=sub2ind([size(mymask,1),size(mymask,2)],globalmerfishspots(:,end-3)-47990+1,globalmerfishspots(:,end-2));
globalmerfishspots(:,13)=ids;
cellstats=regionprops(mymask,'Area','PixelIdxList');
for i=1:size(cellstats,1)
    mypixels=cellstats(i).PixelIdxList;
    positions=ismember(globalmerfishspots(:,13),mypixels);
    globalmerfishspots(positions,14)=i;%add a cell ID to each transcript
end
countable=zeros(size(cellstats,1),3+max(globalmerfishspots(:,2)));
for thiscell=1:size(cellstats,1)%build the count table
    if cellstats(thiscell).Area>0 
        countable(thiscell,1)=thiscell;%#index of that cell
        countable(thiscell,2)=cellstats(thiscell).Area;
        mytable2=globalmerfishspots(globalmerfishspots(:,14)==thiscell,:);
        for gene=1:max(globalmerfishspots(:,2))+1
            cellcountable(1,gene)=sum(mytable2(:,2)==gene-1);
        end
        countable(thiscell,3:end)=cellcountable;
    end
end
writematrix(countable,fullfile(matlabtempfiles,'merfishcounttableregion1rightindex.csv'))
