matlabtempfiles='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/matlabtempfiles/'
%% region 0
calciumpath='/broad/clearylab/Users/Loic/astrocytesFeb2/dish2/20240202_104358_798well2'

matlabmerfish='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/matlabtempfiles'
transformpath='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/outputs'
load(fullfile(matlabmerfish,'mymaskregion0V2.mat'))
mymask=bwlabel(mymask>0);
mymask=downsample(mymask',30);
mymask=downsample(mymask',30);
mymask=flipdim(mymask,1);
load(fullfile(transformpath,'tformMerfishstraighttocalciumregion0well2.mat'))
mymask=imrotate(mymask,-90);
mymask = imwarp_same(mymask,tform);
imwrite(mymask,fullfile(matlabtempfiles,'preppedmaskR0.tif'))
% recoveredmerfish=recoveredmerfish(1:1480,1:1480);%calcium images are 1480x1480
% recoveredmerfish=recoveredmerfish(:,3100:5000);
% recoveredmerfish=[zeros(size(recoveredmerfish,1),193),recoveredmerfish];
% recoveredmerfish=recoveredmerfish(:,3100:5000);
% imwrite(recoveredmerfish,fullfile(matlabtempfiles,'mappedmerfish.tif'))

numberofcells=max(max(mymask));

%loadstack
numberofimages=249
calciumsignal=zeros(numberofcells,numberofimages);

% load(fullfile(matlabtempfiles,'tformMtoCfromalginedshrunkandroratedMregion0.mat'))
for calciumimage=1:numberofimages
    if calciumimage==1
        myimage=imrotate(imread(fullfile(calciumpath,'Time00000_ChannelSS_FITC1_Seq0000.tiff')),180);
        mymask=mymask(1:size(myimage,1),1:size(myimage,2));
    else
        if calciumimage<10
            myimage=imrotate(imread(fullfile(calciumpath,strcat('Time0000',num2str(calciumimage),'_ChannelSS_FITC1_Seq000',num2str(calciumimage),'.tiff'))),180);
        elseif calciumimage<100
            myimage=imrotate(imread(fullfile(calciumpath,strcat('Time000',num2str(calciumimage),'_ChannelSS_FITC1_Seq00',num2str(calciumimage),'.tiff'))),180);
        elseif calciumimage<1000
            myimage=imrotate(imread(fullfile(calciumpath,strcat('Time00',num2str(calciumimage),'_ChannelSS_FITC1_Seq0',num2str(calciumimage),'.tiff'))),180);
        end
    end
    for thiscell=1:numberofcells
        pixelIds=find(mymask==thiscell);
        if max(pixelIds)<size(myimage,1)*size(myimage,2)
            calciumsignal(thiscell,calciumimage)=mean(myimage(pixelIds));
        end
    end
end
disp('finished 1st table')
imwrite(mymask,fullfile(matlabtempfiles,'croppedToCalciumMaskR0.tif'))
imwrite(myimage/100,fullfile(matlabtempfiles,'mappedcalciumR0well2.tif'))
writematrix(calciumsignal,fullfile(matlabtempfiles,'calciumsignalASD2region0savedmemory.csv'))
%% region 1
matlabtempfiles='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/matlabtempfiles/'
calciumpath='/broad/clearylab/Users/Loic/astrocytesFeb2/dish2/20240202_104039_773well1'

matlabmerfish='/broad/clearylab/Users/Loic/merfishdish2/202402142122_loicasdRevsDish2_VMSC11302/matlabtempfiles'
transformpath='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/outputs'
load(fullfile(matlabmerfish,'mymaskregion1V2.mat'))
mymask=bwlabel(mymask>0);
mymask=downsample(mymask',30);
mymask=downsample(mymask',30);
mymask=flipdim(mymask,1);
load(fullfile(transformpath,'tformMerfishstraighttocalciumregion1well1.mat'))
mymask=imrotate(mymask,-90);
mymask = imwarp_same(mymask,tform);
imwrite(mymask,fullfile(matlabtempfiles,'preppedmaskR1.tif'))
% recoveredmerfish=recoveredmerfish(1:1480,1:1480);%calcium images are 1480x1480
% recoveredmerfish=recoveredmerfish(:,750:2500);

numberofcells=max(max(mymask));

%loadstack
numberofimages=249
calciumsignal=zeros(numberofcells,numberofimages);
for calciumimage=1:numberofimages
    if calciumimage==1
        myimage=imrotate(imread(fullfile(calciumpath,'Time00000_ChannelSS_FITC1_Seq0000.tiff')),180);
        mymask=mymask(1:size(myimage,1),1:size(myimage,2));
    else
        if calciumimage<10
            myimage=imrotate(imread(fullfile(calciumpath,strcat('Time0000',num2str(calciumimage),'_ChannelSS_FITC1_Seq000',num2str(calciumimage),'.tiff'))),180);
        elseif calciumimage<100
            myimage=imrotate(imread(fullfile(calciumpath,strcat('Time000',num2str(calciumimage),'_ChannelSS_FITC1_Seq00',num2str(calciumimage),'.tiff'))),180);
        elseif calciumimage<1000
            myimage=imrotate(imread(fullfile(calciumpath,strcat('Time00',num2str(calciumimage),'_ChannelSS_FITC1_Seq0',num2str(calciumimage),'.tiff'))),180);
        end
    end
    for thiscell=1:numberofcells
        pixelIds=find(mymask==thiscell);
        if max(pixelIds)<size(myimage,1)*size(myimage,2)
            calciumsignal(thiscell,calciumimage)=mean(myimage(pixelIds));
        end
    end
end
disp('finished 2nd table')
imwrite(mymask,fullfile(matlabtempfiles,'croppedToCalciumMaskR1.tif'))
imwrite(sum(myimage(:,:,1:30)/1000,3),fullfile(matlabtempfiles,'mappedcalciumR1well1.tif'))
writematrix(calciumsignal,fullfile(matlabtempfiles,'calciumsignalASD2region1savingmemory.csv'))