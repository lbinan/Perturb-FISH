matlabtempfiles='\\helium\broad_clearylab\Users\Loic\tumorMERFISHdish1Donor1FFPE\202408041743_tumor1DonorMerfishFFPE_VMSC11302\matlabtempfiles'
merfishObject = matfile(fullfile(matlabtempfiles,'registrationimageMERFISH.mat'));
merfish1 = double(merfishObject.registrationimage(120000:169973,1:30000)); 
merfish1=reduceImage(reduceImage(merfish1));

merfish1=flipdim(merfish1,1);
figure, imshow(merfish1)

merfish2 = double(merfishObject.registrationimage(120000:169973,30001:50000)); 
merfish2=reduceImage(reduceImage(merfish2));

merfish2=flipdim(merfish2,1);
merfish=[merfish1,merfish2];
figure, imshow(merfish)

% maskObject= matfile(fullfile(matlabtempfiles,'mymaskregion0V2.mat'));
% mymask = double(maskObject.mymask(1:51001,1:10000)); %left part of S2R1
% mymask=flipdim(mymask,1);
% figure, imshow(mymask)


matlabtempfiles='\\helium\broad_clearylab\Users\Loic\tumorZombiedish1Donor1FFPE\matlabtempfiles'
zombieObject = matfile(fullfile(matlabtempfiles,'registrationimageZOMBIE.mat'));
zombie1 = double(zombieObject.registrationimage(1:40000,1:20000)); %for S2R1
zombie1=reduceImage(reduceImage(zombie1));
zombie2 = double(zombieObject.registrationimage(1:40000,20000:50000)); %for S2R1
zombie2=reduceImage(reduceImage(zombie2));
zombie=[zombie1,zombie2];

figure, imshow(zombie)

[movingPoints,fixedPoints] = cpselect(merfish,zombie,"Wait",true);
tform = fitgeotrans(movingPoints,fixedPoints,"similarity");
tformother = fitgeotrans(movingPoints,fixedPoints,"nonreflectivesimilarity");
recovered = imwarp_same(merfish,tform);
% recoveredTranlated=imtranslate(recovered,[-150,20]);
figure, imshowpair(imadjust(zombie),imadjust(recovered))
save(fullfile(matlabtempfiles,'tranformMtoZcompressed.mat'),'tform')

%% test x4
tformlarge=tform;
% save(fullfile(matlabtempfiles,'tranformMtoZReal.mat'),'tformlarge')
merfish1 = double(merfishObject.registrationimage(130000:169973,1:30000)); %bottom left part of S2R1, becomes top left after flip
merfish1=flipdim(merfish1,1);
zombie1 = double(zombieObject.registrationimage(1:15000,1:15000)); %for S2R1
recovered = imwarp_same(merfish1,tformlarge);
% figure, imshowpair(zombie1(1:30000,1:10000),imadjust(recovered))
figure, imshowpair(imadjust(zombie1),imadjust(recovered))

tform=tformlarge;
save(fullfile(matlabtempfiles,'tranformMtoZReal.mat'),'tform','-v7.3')

%% region 1
matlabtempfiles='\\helium\broad_clearylab\Users\Loic\merfishdish2\202402142122_loicasdRevsDish2_VMSC11302\matlabtempfiles'
merfishObject = matfile(fullfile(matlabtempfiles,'cellsmosaicregion1.mat'));
merfish1 = double(merfishObject.cellsmosaic(45000:85600,1:10000)); %bottom left part of S2R1, becomes top left after flip
merfish1=reduceImage(reduceImage(merfish1));
merfish2 = double(merfishObject.cellsmosaic(45000:85600,10000:20000)); %bottom left part of S2R1, becomes top left after flip
merfish2=reduceImage(reduceImage(merfish2));
merfish=[merfish1,merfish2];
merfish=flipdim(merfish,1);
figure, imshow(merfish)

% maskObject= matfile(fullfile(matlabtempfiles,'mymaskregion0V2.mat'));
% mymask = double(maskObject.mymask(1:51001,1:10000)); %left part of S2R1
% mymask=flipdim(mymask,1);
% figure, imshow(mymask)


matlabtempfiles='\\helium\broad_clearylab\Users\Loic\loicAsdRev2_zombie\matlabtempfiles'
zombieObject = matfile(fullfile(matlabtempfiles,'cellsmosaiczombie.mat'));
zombie1 = double(zombieObject.cellsmosaic(60000:80000,1:20000)); %for S2R1
zombie1=reduceImage(reduceImage(zombie1));
zombie2 = double(zombieObject.cellsmosaic(80000:100000,1:20000)); %for S2R1
zombie2=reduceImage(reduceImage(zombie2));
zombie=[zombie1;zombie2];

figure, imshow(zombie)

[movingPoints,fixedPoints] = cpselect(merfish,zombie/3,"Wait",true);
tform = fitgeotrans(movingPoints,fixedPoints,"similarity");
% tformother = fitgeotrans(movingPoints,fixedPoints,"nonreflectivesimilarity");
recovered = imwarp_same(merfish,tform);
% recoveredTranlated=imtranslate(recovered,[-150,20]);
figure, imshowpair(imadjust(zombie),imadjust(recovered))

%% test x4
tformlarge=tform;
merfish1 = double(merfishObject.cellsmosaic(100000:151000,1:10000)); %bottom left part of S2R1, becomes top left after flip
merfish1=flipdim(merfish1,1);
zombie1 = double(zombieObject.cellsmosaic(1:30000,1:10000)); %for S2R1
recovered = imwarp_same(merfish1(1:30000,1:10000),tformlarge);
% figure, imshowpair(zombie1(1:30000,1:10000),imadjust(recovered))
figure, imshowpair(imadjust(zombie1(1:30000,1:10000)),imadjust(recovered))

tform=tformlarge;
save(fullfile(matlabtempfiles,'tformMerfS2R1TOzombietopFromreducedImages.mat'),'tform','-v7.3')



% load(fullfile(matlabtempfiles,'tformMerfS2R1TOzombietopFromreducedImages.mat'))



%%
%% sample4
%%

matlabtempfiles='\\helium\broad_clearylab\Users\Loic\ASDMerfishRevs4\matlabtempfiles'
merfishObject = matfile(fullfile(matlabtempfiles,'cellsmosaicregion0.mat'));
merfish1 = double(merfishObject.cellsmosaic(64322:133760,1:10000)); %bottom left part of S2R1, becomes top left after flip
merfish1=reduceImage(reduceImage(merfish1));
merfish2 = double(merfishObject.cellsmosaic(64322:133760,10000:20000)); %bottom left part of S2R1, becomes top left after flip
merfish2=reduceImage(reduceImage(merfish2));
merfish=[merfish1,merfish2];
merfish=flipdim(merfish,1);
figure, imshow(merfish)

% maskObject= matfile(fullfile(matlabtempfiles,'mymaskregion0V2.mat'));
% mymask = double(maskObject.mymask(1:51001,1:10000)); %left part of S2R1
% mymask=flipdim(mymask,1);
% figure, imshow(mymask)


matlabtempfiles='\\helium\broad_clearylab\Users\Loic\ASdRevsZombiedish4\matlabtempfiles'
zombieObject = matfile(fullfile(matlabtempfiles,'cellsmosaiczombie.mat'));
zombie1 = double(zombieObject.cellsmosaic(1:20000,1:20000)); %for S2R1
zombie1=reduceImage(reduceImage(zombie1));
zombie2 = double(zombieObject.cellsmosaic(20000:50000,1:20000)); %for S2R1
zombie2=reduceImage(reduceImage(zombie2));
zombie=[zombie1;zombie2];

figure, imshow(zombie)

[movingPoints,fixedPoints] = cpselect(merfish,zombie/3,"Wait",true);
tform = fitgeotrans(movingPoints,fixedPoints,"similarity");
tformother = fitgeotrans(movingPoints,fixedPoints,"nonreflectivesimilarity");
recovered = imwarp_same(merfish,tform);
% recoveredTranlated=imtranslate(recovered,[-150,20]);
figure, imshowpair(imadjust(zombie),imadjust(recovered))

%% test x4
tformlarge=tform;
merfish1 = double(merfishObject.cellsmosaic(84322:133760,1:10000)); %bottom left part of S2R1, becomes top left after flip
merfish1=flipdim(merfish1,1);
zombie1 = double(zombieObject.cellsmosaic(1:30000,1:10000)); %for S2R1
recovered = imwarp_same(merfish1(1:30000,1:10000),tformlarge);
% figure, imshowpair(zombie1(1:30000,1:10000),imadjust(recovered))
figure, imshowpair(imadjust(zombie1(1:30000,1:10000)),imadjust(recovered))

tform=tformlarge;
save(fullfile(matlabtempfiles,'tformMerfS4TOzombietopFromreducedImages.mat'),'tform','-v7.3')
