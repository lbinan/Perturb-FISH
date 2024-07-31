% % % mypath='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/data'
% % % matlabtempfiles='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/matlabtempfiles'
% % % nucleipath='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/nuclei'
% % % %% region0
% % % dapiObject = matfile(fullfile(matlabtempfiles,'binaryDapiregion0.mat'));
% % % % size(dapiObject.dapimosaic) real image is 167925 by 75765
% % % % dapiObject.Properties.Writable = true;
% % % dapi = dapiObject.dapimosaic(98800:151000,:); %for S2R0
% % % dapi=dapi>0;
% % % cells=regionprops(dapi,'Centroid');
% % % seeds=zeros(size(dapi,1),size(dapi,2));
% % % for i=1:size(cells,1)
% % %     seeds(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
% % % end
% % % seeds=seeds>0;
% % %     
% % % largedapi=imdilate(dapi, strel('disk',50));
% % % D = bwdist(largedapi);
% % % gmag2 = imimposemin(D , seeds);
% % % mymask = watershed(gmag2);
% % % mymask(~largedapi)=0;
% % % mymask=mymask>0;
% % % save(fullfile(matlabtempfiles,'mymaskregion0V2.mat'),'mymask','-v7.3')
% % % clear
% % % mypath='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/data'
% % % matlabtempfiles='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/matlabtempfiles'
% % % nucleipath='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/nuclei'
% % % %% region1
% % % dapiObject = matfile(fullfile(matlabtempfiles,'binaryDapiregion1.mat'));
% % % % size(dapiObject.dapimosaic) real image is 94208 by 75755
% % % % dapiObject.Properties.Writable = true;
% % % dapi = double(dapiObject.dapimosaic(32000:85600,1:67200)); %for S2R1
% % % dapi=dapi>0;
% % % cells=regionprops(dapi,'Centroid');
% % % seeds=zeros(size(dapi,1),size(dapi,2));
% % % for i=1:size(cells,1)
% % %     seeds(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
% % % end
% % % seeds=seeds>0;
% % %     
% % % largedapi=imdilate(dapi, strel('disk',50));
% % % D = bwdist(largedapi);
% % % gmag2 = imimposemin(D , seeds);
% % % mymask = watershed(gmag2);
% % % mymask(~largedapi)=0;
% % % mymask=mymask>0;
% % % 
% % % save(fullfile(matlabtempfiles,'mymaskregion1V2.mat'),'mymask','-v7.3')
% % % clear
% % % 
% % % %%
% % % disp('next_sample')
% % % %%
% % % mypath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/data'
% % % matlabtempfiles='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/matlabtempfiles'
% % % nucleipath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/nuclei'
% % % %% region0
% % % dapiObject = matfile(fullfile(matlabtempfiles,'binaryDapiregion0.mat'));
% % % % size(dapiObject.dapimosaic) real image is 122880 by 77813
% % % % dapiObject.Properties.Writable = true;
% % % dapi = double(dapiObject.dapimosaic(41600:112000,12206:71190)); %for S3R0
% % % dapi=dapi>0;
% % % cells=regionprops(dapi,'Centroid');
% % % seeds=zeros(size(dapi,1),size(dapi,2));
% % % for i=1:size(cells,1)
% % %     seeds(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
% % % end
% % % seeds=seeds>0;
% % %     
% % % largedapi=imdilate(dapi, strel('disk',50));
% % % D = bwdist(largedapi);
% % % gmag2 = imimposemin(D , seeds);
% % % mymask = watershed(gmag2);
% % % mymask(~largedapi)=0;
% % % mymask=mymask>0;
% % % save(fullfile(matlabtempfiles,'mymaskregion0V2.mat'),'mymask','-v7.3')
% % % clear
% % % mypath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/data'
% % % matlabtempfiles='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/matlabtempfiles'
% % % nucleipath='/broad/clearylab/Users/Loic/asdrevsdish3merfish/202402171543_loicAsdRevsDish3_VMSC11302/nuclei'
% % % %% region1
% % % dapiObject = matfile(fullfile(matlabtempfiles,'binaryDapiregion1.mat'));
% % % % size(dapiObject.dapimosaic) real image is 129013 by 63477
% % % % dapiObject.Properties.Writable = true;
% % % dapi = double(dapiObject.dapimosaic(47990:116000,1:54400)); %for S3R1
% % % dapi=dapi>0;
% % % load(fullfile(matlabtempfiles,'globalmerfishspotsregion1.mat'));
% % % cells=regionprops(dapi,'Centroid');
% % % seeds=zeros(size(dapi,1),size(dapi,2));
% % % for i=1:size(cells,1)
% % %     seeds(floor(cells(i).Centroid(2)),floor(cells(i).Centroid(1)))=1;
% % % end
% % % seeds=seeds>0;
% % %     
% % % largedapi=imdilate(dapi, strel('disk',50));
% % % D = bwdist(largedapi);
% % % gmag2 = imimposemin(D , seeds);
% % % mymask = watershed(gmag2);
% % % mymask(~largedapi)=0;
% % % mymask=mymask>0;
% % % save(fullfile(matlabtempfiles,'mymaskregion1V2.mat'),'mymask','-v7.3')

%%
clear
revisionsmerfishcounttable2
clear
revisionsmerfishcounttable3


