mypath='/broad/clearylab/Users/Loic/ASDMerfishRevs4/analysed/data'
matlabtempfiles='/broad/clearylab/Users/Loic/ASDMerfishRevs4/matlabtempfiles'
nucleipath='/broad/clearylab/Users/Loic/ASDMerfishRevs4/nuclei'
load(fullfile(matlabtempfiles,'mymaskregion0V2.mat'));
load(fullfile(matlabtempfiles,'globalmerfishspotsregion0.mat'));
globalmerfishspots=globalmerfishspots(globalmerfishspots(:,end-3)<133760,:);
% for i=1:size(globalmerfishspots,1)
%     spotsimage(floor(globalmerfishspots(i,end-3)),floor(globalmerfishspots(i,end-2)))=globalmerfishspots(i,1);
% end
mymask=mymask>0;
%  ids=find(spotsimage)>0;
ids=sub2ind([size(mymask,1),size(mymask,2)],globalmerfishspots(:,end-3)-64322+1,globalmerfishspots(:,end-2)+1);
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
writematrix(countable,fullfile(matlabtempfiles,'merfishcounttablerightIndex.csv'))
%%