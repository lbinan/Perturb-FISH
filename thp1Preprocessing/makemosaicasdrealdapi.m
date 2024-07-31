localpath='\\helium\broad_clearylab\Users\Loic\thp1homemadezombie_1\Data'
globalpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/Data'
savepath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/preprocessed/'
matlabtempfiles='/broad/clearylab/Users/Loic/thp1homemadezombie_1/matlabtempfiles/'
maxprojpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/maxproj/'
filteredpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/filteredimages/'
localpath=globalpath;
minX=0;
minY=0;
maxX=0;
maxY=0;
myfovs=[0:1275];
for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    if thisFOV<10
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_1.xml')));
    elseif thisFOV<100
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_1.xml')));
    elseif thisFOV<1000
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_1.xml')));
    else
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_1.xml')));
    end
    thisposition=thisinfo.settings.acquisition.stage_position.Text;
    C = textscan(thisposition, '%f,%f');
    thisX=C{1}/200;
    thisY=C{2}/200;
    minX=min(minX,thisX);%finding size of mosaic in number of images, and finding the center
    minY=min(minY,thisY);
    maxX=max(maxX,thisX);
    maxY=max(maxY,thisY);
end
% mymosaic=imread(fullfile(matlabtempfiles,'testmosaic4.dax'));
milieuX=round((abs(minX)+abs(maxX))/2)+1;%finding the center
milieuY=round((abs(minY)+abs(maxY))/2)+1;

mymosaic=im2uint8(zeros(uint16(1024+(maxX-minX)*1024),uint16(1024+(maxY-minY)*1024)));
for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    if thisFOV<10
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_000',num2str(thisFOV),'_1.xml')));
    %     thisname=fullfile(maxprojpath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_000',num2str(thisFOV),'.dax'));
        myimage=ReadDax(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_000',num2str(thisFOV),'.dax')));
        myimage = im2double(myimage(:,:,30));%loading only dapi
    elseif thisFOV<100
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_00',num2str(thisFOV),'_1.xml')));
    % 	thisname=fullfile(maxprojpath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_00',num2str(thisFOV),'.dax'));
        myimage=ReadDax(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_00',num2str(thisFOV),'.dax')));
        myimage = im2double(myimage(:,:,30));
    elseif thisFOV<1000
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_0',num2str(thisFOV),'_1.xml')));
    %     thisname=fullfile(maxprojpath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_0',num2str(thisFOV),'.dax'));
        myimage=ReadDax(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_0',num2str(thisFOV),'.dax')));
        myimage = im2double(myimage(:,:,30));
    else
        thisinfo=xml2struct(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z1-405z1_',num2str(thisFOV),'_1.xml')));
        myimage=ReadDax(fullfile(localpath,strcat('hal-config-749z7-638z7-546z7-477z7-405z7_',num2str(thisFOV),'.dax')));
        myimage = im2double(myimage(:,:,30));
    end
    % thisimage=im2uint8(max(myimage(:,:,1:7),[],3));
    thisimage=im2uint8(reduceImage(myimage(:,:,1)));%rescale, and correct the orientation of the camera on the scope
    thisimage=flipdim(thisimage,1);
    thisimage=flipdim(thisimage,2);
    thisimage=imrotate(thisimage,-90);
    thisposition=thisinfo.settings.acquisition.stage_position.Text;
    C = textscan(thisposition, '%f,%f');
    thisX=C{1}/200+abs(minX)+1;
    thisY=C{2}/200+abs(minY)+1;
    %insert image into mosaic. taking max deals with overlapping region
    mymosaic((thisY-1)*922+1:(thisY-1)*922+1024,(thisX-1)*922+1:(thisX-1)*922+1024)=max(mymosaic((thisY-1)*922+1:(thisY-1)*922+1024,(thisX-1)*922+1:(thisX-1)*922+1024),thisimage);
end
mymosaic=mymosaic+1;
imwrite(mymosaic,fullfile(matlabtempfiles,strcat('testrealdapi.tif')))
save(fullfile(matlabtempfiles,'dapimosaic.mat'),'mymosaic','-v7.3')