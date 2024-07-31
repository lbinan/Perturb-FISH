matlabtempfiles='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/matlabtempfiles'
%% region0
dapiObject = matfile(fullfile(matlabtempfiles,'binaryDapiregion0.mat'));
% size(dapiObject.dapimosaic) real image is 167925 by 75765
% dapiObject.Properties.Writable = true;
dapi = dapiObject.dapimosaic(41600:112000,12206:71190); %for S3R0
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
matlabtempfiles='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/matlabtempfiles'
%% region1
dapiObject = matfile(fullfile(matlabtempfiles,'binaryDapiregion1.mat'));
% size(dapiObject.dapimosaic) real image is 94208 by 75755
% dapiObject.Properties.Writable = true;
dapi = double(dapiObject.dapimosaic(47990:116000,1:54400)); %for S3R1
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

save(fullfile(matlabtempfiles,'mymaskregion1V2.mat'),'mymask','-v7.3')
clear

%%