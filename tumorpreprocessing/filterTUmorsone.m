mypath='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/Data'

filtered='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/filtered2'
matlabtempfiles='/broad/clearylab/Users/Loic/tumorZombiedish1Donor1FFPE/matlabtempfiles'
blankpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1'

localpath=mypath
myfovs=[0:1169];
delete(gcp('nocreate'))

thisfilesname='blank_0001.dax'
background1=double(ReadDax(fullfile(blankpath,thisfilesname)));
thisfilesname='test_0002.dax'
background2=double(ReadDax(fullfile(blankpath,thisfilesname)));
for i=1:21
    background(:,:,i)=(background1(:,:,i)+background2(:,:,i))/2;
end
clear background1 background2


% myparpool=parpool(22)
% parfor (fov=1:size(myfovs,2),22)
for fov=1:size(myfovs,2)
    filteredimage=zeros(2048,2048,21);%rebuild file name, for cycle -1 
    %(the scope acquire images in cycle, and counts from -1, then 0 then 1 etc...
    %images for round :-1" have slighly different names because they are a
    %larger stack stat countains more images for DAPI and registration
    %beads
    thisFOV=myfovs(fov);
    if thisFOV<10
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z7-405z7_000',num2str(thisFOV),'.dax');
    elseif thisFOV<100
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z7-405z7_00',num2str(thisFOV),'.dax');
    elseif thisFOV<1000
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z7-405z7_0',num2str(thisFOV),'.dax');
    else
        thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z7-405z7_',num2str(thisFOV),'.dax');
    end
    imageref=double(ReadDax(fullfile(localpath,thisfilesname)));
    for thiscycle=1:5
         if thisFOV<10
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');
        elseif thisFOV<100
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');
        elseif thisFOV<1000
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');
        else
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_',num2str(thiscycle-1),'.dax');
         end
         if thiscycle==1
        if thisFOV<10
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_',num2str(thiscycle),'.dax');
        elseif thisFOV<100
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_',num2str(thiscycle),'.dax');
        elseif thisFOV<1000
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_',num2str(thiscycle),'.dax');
        else
            thisfilesname=strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_',num2str(thiscycle),'.dax');
        end
         end
        thatfilename=thisfilesname;
            myimage=ReadDax(fullfile(localpath,thisfilesname));
        if thiscycle==1
            imageround1=registerImages(imageref,double(myimage));
        elseif thiscycle==2
            imageround2=registerImages(imageref,double(myimage));
        elseif thiscycle==3
            imageround3=registerImages(imageref,double(myimage));
        elseif thiscycle==4
            for pre=1:size(myimage,3)
                myimage(:,:,pre)=imtranslate(double(myimage(:,:,pre)),[-21,96]);
            end
            imageround4=registerImages(imageref,myimage);
        elseif thiscycle==5
            for pre=1:size(myimage,3)
                myimage(:,:,pre)=imtranslate(double(myimage(:,:,pre)),[-21,96]);
            end
            imageround5=registerImages(imageref,myimage);
        end
    end
    for i=1:21%remove non uniform illumination by dividing by background
        imageround1(:,:,i)=imageround1(:,:,i)./background(:,:,i);
        imageround2(:,:,i)=imageround2(:,:,i)./background(:,:,i);
        imageround3(:,:,i)=imageround3(:,:,i)./background(:,:,i);
        imageround4(:,:,i)=imageround4(:,:,i)./background(:,:,i);
        imageround5(:,:,i)=imageround5(:,:,i)./background(:,:,i);
    end
    for frame=1:21%gaussian filter to enhance spots of the expected size for a guide
        imageround1(:,:,frame)=imgaussfilt(imageround1(:,:,frame),2)-imgaussfilt(imageround1(:,:,frame),12);
        imageround2(:,:,frame)=imgaussfilt(imageround2(:,:,frame),2)-imgaussfilt(imageround2(:,:,frame),12);
        imageround3(:,:,frame)=imgaussfilt(imageround3(:,:,frame),2)-imgaussfilt(imageround3(:,:,frame),12);
        imageround4(:,:,frame)=imgaussfilt(imageround4(:,:,frame),2)-imgaussfilt(imageround4(:,:,frame),12);
        imageround5(:,:,frame)=imgaussfilt(imageround5(:,:,frame),2)-imgaussfilt(imageround5(:,:,frame),12);
    end
    croppedfilename=thisfilesname(1:end-5);%make a new file name
    %save max projected images, projecting on Z for images of a given
    %channel(=color)
    imwrite(max(imageround1(:,:,2:7),[],3)/2,fullfile(filtered,strcat(croppedfilename,'1.tif')));
    imwrite(max(imageround1(:,:,8:13),[],3)/2,fullfile(filtered,strcat(croppedfilename,'2.tif')));
    imwrite(max(imageround1(:,:,16:21),[],3)/2,fullfile(filtered,strcat(croppedfilename,'3.tif')));
    imwrite(max(imageround2(:,:,2:7),[],3)/2,fullfile(filtered,strcat(croppedfilename,'4.tif')));
    imwrite(max(imageround2(:,:,8:13),[],3)/2,fullfile(filtered,strcat(croppedfilename,'5.tif')));
    imwrite(max(imageround2(:,:,16:21),[],3)/2,fullfile(filtered,strcat(croppedfilename,'6.tif')));
    imwrite(max(imageround3(:,:,2:7),[],3)/2,fullfile(filtered,strcat(croppedfilename,'7.tif')));
    imwrite(max(imageround3(:,:,8:13),[],3)/2,fullfile(filtered,strcat(croppedfilename,'8.tif')));
    imwrite(max(imageround3(:,:,16:21),[],3)/2,fullfile(filtered,strcat(croppedfilename,'9.tif')));
    imwrite(max(imageround4(:,:,2:7),[],3)/2,fullfile(filtered,strcat(croppedfilename,'10.tif')));
    imwrite(max(imageround4(:,:,8:13),[],3)/2,fullfile(filtered,strcat(croppedfilename,'11.tif')));
    imwrite(max(imageround4(:,:,16:21),[],3)/2,fullfile(filtered,strcat(croppedfilename,'12.tif')));
    imwrite(max(imageround5(:,:,2:7),[],3)/2,fullfile(filtered,strcat(croppedfilename,'13.tif')));
    imwrite(max(imageround5(:,:,8:13),[],3)/2,fullfile(filtered,strcat(croppedfilename,'14.tif')));
    imwrite(max(imageround5(:,:,16:21),[],3)/4,fullfile(filtered,strcat(croppedfilename,'15.tif')));
end
    
delete(gcp('nocreate'))

disp(('gothere'))
decodetumors

%%
% disp('start mosaic')
% maskeZOMBIEmosaictumor1
% disp('did mosaic')
