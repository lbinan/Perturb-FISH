mypath='/broad/clearylab/Users/Loic/ASDMerfishRevs4/analysed/data'
matlabtempfiles='/broad/clearylab/Users/Loic/ASDMerfishRevs4/matlabtempfiles'
nucleipath='/broad/clearylab/Users/Loic/ASDMerfishRevs4/nuclei'
%% region0
dapiObject = matfile(fullfile(matlabtempfiles,'binaryDapiregion0.mat'));
% size(dapiObject.dapimosaic) real image is 147445 by 75765
% dapiObject.Properties.Writable = true;
dapi = double(dapiObject.dapimosaic(64322:133760,:)); %for S4R0
dapi=dapi>0;
cells=regionprops(dapi,'Centroid');
seeds=zeros(size(dapi,1),size(dapi,2));
for i=1:size(cells,1)
    seeds(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
end
seeds=seeds>0;
    
largedapi=imdilate(dapi, strel('disk',50));
D = bwdist(largedapi);
gmag2 = imimposemin(D , seeds);
mymask = watershed(gmag2);
mymask(~largedapi)=0;
mymask=mymask>0;
save(fullfile(matlabtempfiles,'mymaskregion0V2.mat'),'mymask','-v7.3')
clear

% % clear
% % revisionsmerfishcounttable2
% % clear
% % revisionsmerfishcounttable3
% % 
% revisionsmerfishcounttable4
