mypath='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/Data'%path to raw images (dax files)
filtered='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/filtered'%output path
matlabtempfiles='/broad/clearylab/Users/Loic/loicAsdRev2_zombie/matlabtempfiles'
blankpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1'%path to 2 images of a blank area, to compensate for non-uniform illumination

localpath=mypath
myfovs=[0:2010];
delete(gcp('nocreate'))

thisfilesname='blank_0001.dax'
background1=double(ReadDax(fullfile(blankpath,thisfilesname)));
thisfilesname='test_0002.dax'
background2=double(ReadDax(fullfile(blankpath,thisfilesname)));
for i=1:21
    background(:,:,i)=(background1(:,:,i)+background2(:,:,i))/2;%average the 2 blank image stacks
end
clear background1 background2


myparpool=parpool(22)
parfor (fov=1:size(myfovs,2),22)
% filteredimage=zeros(2048,2048,16);
% imageround1=zeros(2048,2048,16);
% imageround2=zeros(2048,2048,16);
% imageround3=zeros(2048,2048,16);
% imageround4=zeros(2048,2048,16);
% imageround5=zeros(2048,2048,16);
% imageround6=zeros(2048,2048,16);
% imageround7=zeros(2048,2048,16);
% imageround8=zeros(2048,2048,16);
% imageround9=zeros(2048,2048,16);
% imageround10=zeros(2048,2048,16);
% imageround11=zeros(2048,2048,16);
thisFOV=myfovs(fov);
if thisFOV<10
    thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z7-405z7_000',num2str(thisFOV),'.dax');
    otherfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_000',num2str(thisFOV),'_00.dax');
elseif thisFOV<100
    thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z7-405z7_00',num2str(thisFOV),'.dax');
    otherfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_00',num2str(thisFOV),'_00.dax');
elseif thisFOV<1000
    thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z7-405z7_0',num2str(thisFOV),'.dax');
    otherfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_0',num2str(thisFOV),'_00.dax');
else
    thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z7-405z7_',num2str(thisFOV),'.dax');
    otherfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_',num2str(thisFOV),'_00.dax');
end
croppedfilename=otherfilesname(1:end-6);%rename
% if ~isfile(fullfile(filtered,strcat(croppedfilename,'1.tif'))) | ~isfile(fullfile(filtered,strcat(croppedfilename,'21.tif')))

imageround0=double(ReadDax(fullfile(localpath,thisfilesname)));
for thiscycle=1:11
     if thisFOV<10
         if thiscycle<11
             thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_000',num2str(thisFOV),'_0',num2str(thiscycle-1),'.dax');
         else
             thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_000',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');
         end
    elseif thisFOV<100
         if thiscycle<11
             thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_00',num2str(thisFOV),'_0',num2str(thiscycle-1),'.dax');
         else
             thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_00',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');
         end
    elseif thisFOV<1000
         if thiscycle<11
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_0',num2str(thisFOV),'_0',num2str(thiscycle-1),'.dax');
         else
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_0',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');
         end
    else
         if thiscycle<11
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_',num2str(thisFOV),'_0',num2str(thiscycle-1),'.dax');
         else
            thisfilesname=strcat('bicolorhal-config-749z7-638z7-477z1-405z1_',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');   
         end
     end
    thatfilename=thisfilesname;
            
    myimage=ReadDax(fullfile(localpath,thisfilesname));
    if thiscycle==1
        imageround1=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==2
        imageround2=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==3
        imageround3=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==4
        imageround4=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==5
        imageround5=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==6
        imageround6=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==7
        imageround7=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==8
        imageround8=registerRevsImages(imageround0,double(myimage));
    elseif thiscycle==9
        imageround9=registerRevsImages(imageround0,imtranslate(double(myimage),[10,55]));
    elseif thiscycle==10
        imageround10=registerRevsImages(imageround0,imtranslate(double(myimage),[10,55]));
    elseif thiscycle==11
        imageround11=registerRevsImages(imageround0,imtranslate(double(myimage),[10,55]));
    end
end
for i=1:16%images are divided by the average blank to compensate for the shape of the excitation light
imageround1(:,:,i)=imageround1(:,:,i)./background(:,:,i);
imageround2(:,:,i)=imageround2(:,:,i)./background(:,:,i);
imageround3(:,:,i)=imageround3(:,:,i)./background(:,:,i);
imageround4(:,:,i)=imageround4(:,:,i)./background(:,:,i);
imageround5(:,:,i)=imageround5(:,:,i)./background(:,:,i);
imageround6(:,:,i)=imageround6(:,:,i)./background(:,:,i);
imageround7(:,:,i)=imageround7(:,:,i)./background(:,:,i);
imageround8(:,:,i)=imageround8(:,:,i)./background(:,:,i);
imageround9(:,:,i)=imageround9(:,:,i)./background(:,:,i);
imageround10(:,:,i)=imageround10(:,:,i)./background(:,:,i);
imageround11(:,:,i)=imageround11(:,:,i)./background(:,:,i);
end
for frame=1:16%images are gaussian filtered to enhance signal and remove background
imageround1(:,:,frame)=imgaussfilt(imageround1(:,:,frame),2)-imgaussfilt(imageround1(:,:,frame),12);
imageround2(:,:,frame)=imgaussfilt(imageround2(:,:,frame),2)-imgaussfilt(imageround2(:,:,frame),12);
imageround3(:,:,frame)=imgaussfilt(imageround3(:,:,frame),2)-imgaussfilt(imageround3(:,:,frame),12);
imageround4(:,:,frame)=imgaussfilt(imageround4(:,:,frame),2)-imgaussfilt(imageround4(:,:,frame),12);
imageround5(:,:,frame)=imgaussfilt(imageround5(:,:,frame),2)-imgaussfilt(imageround5(:,:,frame),12);
imageround6(:,:,frame)=imgaussfilt(imageround6(:,:,frame),2)-imgaussfilt(imageround6(:,:,frame),12);
imageround7(:,:,frame)=imgaussfilt(imageround7(:,:,frame),2)-imgaussfilt(imageround7(:,:,frame),12);
imageround8(:,:,frame)=imgaussfilt(imageround8(:,:,frame),2)-imgaussfilt(imageround8(:,:,frame),12);
imageround9(:,:,frame)=imgaussfilt(imageround9(:,:,frame),2)-imgaussfilt(imageround9(:,:,frame),12);
imageround10(:,:,frame)=imgaussfilt(imageround10(:,:,frame),2)-imgaussfilt(imageround10(:,:,frame),12);
imageround11(:,:,frame)=imgaussfilt(imageround11(:,:,frame),2)-imgaussfilt(imageround11(:,:,frame),12);
end
croppedfilename=thisfilesname(1:end-6);%rename
%max project all images from one round before saving
imwrite(max(imageround1(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'1.tif')));
imwrite(max(imageround1(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'2.tif')));

imwrite(max(imageround2(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'3.tif')));
imwrite(max(imageround2(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'4.tif')));

imwrite(max(imageround3(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'5.tif')));
imwrite(max(imageround3(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'6.tif')));

imwrite(max(imageround4(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'7.tif')));
imwrite(max(imageround4(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'8.tif')));

imwrite(max(imageround5(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'9.tif')));
imwrite(max(imageround5(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'10.tif')));

imwrite(max(imageround6(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'11.tif')));
imwrite(max(imageround6(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'12.tif')));

imwrite(max(imageround7(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'13.tif')));
imwrite(max(imageround7(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'14.tif')));

imwrite(max(imageround8(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'15.tif')));
imwrite(max(imageround8(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'16.tif')));

imwrite(max(imageround9(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'17.tif')));
imwrite(max(imageround9(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'18.tif')));

imwrite(max(imageround10(:,:,2:7),[],3)/5,fullfile(filtered,strcat(croppedfilename,'19.tif')));
imwrite(max(imageround10(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'20.tif')));

imwrite(max(imageround11(:,:,8:13),[],3)/5,fullfile(filtered,strcat(croppedfilename,'21.tif')));
end

    
delete(gcp('nocreate'))
disp('got here')
decodeRevsastrocytesdish2
spotsQCastrocyteAsd2revs
clear
mypath='/stanley/opticalprofiling_storage/Loic/AsdRevsDish3Zombie/Data'%path to raw images (dax files)
filtered='/broad/clearylab/Users/Loic/loicAsdRev3_zombie/filtered'%output path
matlabtempfiles='/broad/clearylab/Users/Loic/loicAsdRev3_zombie/matlabtempfiles'
blankpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1'%path to 2 images of a blank area, to compensate for non-uniform illumination

localpath=mypath
myfovs=[0:1800];
delete(gcp('nocreate'))
% filterRevAsd3
decodeRevsastrocytesdish3
spotsQCastrocyteAsd3revs
