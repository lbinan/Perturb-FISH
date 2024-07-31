matlabtempfiles='\\helium\broad_clearylab\Users\Loic\merfishdish2\202402142122_loicasdRevsDish2_VMSC11302\matlabtempfiles'
merfishObject = matfile(fullfile(matlabtempfiles,'cellsmosaicregion0.mat'));
merfish1 = double(merfishObject.cellsmosaic(100000:151000,1:30000)); %bottom left part of S2R1, becomes top left after flip
merfish1=reduceImage(reduceImage(merfish1));
merfish2 = double(merfishObject.cellsmosaic(100000:151000,30000:60000)); %bottom left part of S2R1, becomes top left after flip
merfish2=reduceImage(reduceImage(merfish2));
merfish=[merfish1,merfish2];
merfish=flipdim(merfish,1);
figure, imshow(merfish)
merfish1=flipdim(merfish1,1);
matlabzombie='\\helium\broad_clearylab\Users\Loic\loicAsdRev2_zombie\matlabtempfiles'
load(fullfile(matlabzombie,'tformMerfS2R0TOzombietopFromreducedImages.mat'))
recoveredmerfish = imwarp_same(merfish1,tform);
recoveredmerfish=downsample(recoveredmerfish',30);
recoveredmerfish=downsample(recoveredmerfish',30);
recoveredmerfish=imrotate(recoveredmerfish,-90);
imwrite(recoveredmerfish,fullfile(matlabzombie,'merfishpreparedforcalcium.tif'))
clear merfish
calcium=imread('\\helium\broad_clearylab\Users\Loic\astrocytesFeb2\dish2\well2map.tif');
calcium=imrotate(calcium,180);
% merfish=imread(fullfile(matlabzombie,'merfishpreparedforcalcium.tif'));
merfish=recoveredmerfish(:,3100:5000);
[movingPoints,fixedPoints] = cpselect(calcium,merfish,"Wait",true);
tform = fitgeotrans(movingPoints,fixedPoints,"similarity");
recovered = imwarp_same(calcium,tform);
figure, imshowpair(recovered,imadjust(merfish))
save(fullfile(matlabzombie,'tformMtoCfromalginedshrunkandroratedMregion0.mat'),'tform','-v7.3')
%%
matlabtempfiles='\\helium\broad_clearylab\Users\Loic\merfishdish2\202402142122_loicasdRevsDish2_VMSC11302\matlabtempfiles'
merfishObject = matfile(fullfile(matlabtempfiles,'cellsmosaicregion1.mat'));
merfish = double(merfishObject.cellsmosaic(32000:85600,1:67200)); %left part of S2R1
merfish=flipdim(merfish,1);
figure, imshow(merfish)
matlabzombie='\\helium\broad_clearylab\Users\Loic\loicAsdRev2_zombie\matlabtempfiles'
load(fullfile(matlabzombie,'tformMerfS2R1TOzombiebottom.mat'))
recoveredmerfish = imwarp_same(merfish,tform);
recoveredmerfish=downsample(recoveredmerfish',30);
recoveredmerfish=downsample(recoveredmerfish',30);
recoveredmerfish=imrotate(recoveredmerfish,-90);
imwrite(recoveredmerfish,fullfile(matlabzombie,'merfishpreparedforcalciumregion1.tif'))
clear merfish
calcium=imread('\\helium\broad_clearylab\Users\Loic\astrocytesFeb2\dish2\well1map.tif');
calcium=imrotate(calcium,180);
merfish=imread(fullfile(matlabzombie,'merfishpreparedforcalciumregion1.tif'));
[movingPoints,fixedPoints] = cpselect(calcium,merfish,"Wait",true);
tform = fitgeotrans(movingPoints,fixedPoints,"similarity");
recovered = imwarp_same(calcium,tform);
figure, imshowpair(recovered,merfish)
save(fullfile(matlabzombie,'tformMtoCfromalginedshrunkandroratedMregion1.mat'),'tform','-v7.3')

%%
