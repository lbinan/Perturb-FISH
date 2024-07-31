%converts decoded spots into coordinates of global mosaic. helps with
%regions split between images
localpath='\\helium\broad_clearylab\Users\Loic\thp1homemadezombie_1\Data'
globalpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/Data'
savepath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/preprocessed/'
matlabtempfiles='/broad/clearylab/Users/Loic/thp1homemadezombie_1/matlabtempfiles/'
maxprojpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/maxproj/'
filteredpath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/filteredimages/'
merlinpath='/broad/clearylab/Users/Loic/thp1homemademerfish_1/merlin/Data/GenerateMosaic/images'
maskedagain='/broad/clearylab/Users/Loic/thp1homemadezombie_1/maskedagain'
filtered='/broad/clearylab/Users/Loic/thp1homemadezombie_1/filtered'
mypath='/broad/clearylab/Users/Loic/thp1homemadezombie_1/Data'
localpath=globalpath;
minX=0;
minY=0;
maxX=0;
maxY=0;

zombiebarcodes=table2array(readtable(fullfile(matlabtempfiles,'filteredspotsround3SomeMore.csv')));%load decoded spots after QC filtering
globalzombiespots=zombiebarcodes;
clear zombiebarcodes
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
rotatedtemplate=imrotate(thisimage,-90);
myfovs=[0:1275];
thisimage=zeros(2048,2048);
for i=1:2048*2048
    thisimage(i)=i;
end
reverserotatedtemplate=imrotate(thisimage,90);
%%
for fov=1:size(myfovs,2)%first loop to find size of mosaic and center
    thisFOV=myfovs(fov);
    disp(strcat('doing FOV ', num2str(thisFOV)))
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
    minX=min(minX,thisX);
    minY=min(minY,thisY);
    maxX=max(maxX,thisX);
    maxY=max(maxY,thisY);
end
milieuX=round((abs(minX)+abs(maxX))/2)+1;
milieuY=round((abs(minY)+abs(maxY))/2)+1;
mymosaic=im2uint8(zeros(floor(2048+(maxX-minX)*2048),floor(2048+(maxY-minY)*2048)));
%% 
for fov=1:size(myfovs,2)
    thisFOV=myfovs(fov);
    if rem(thisFOV,100)==0
        disp(strcat('current FOV is ',num2str(thisFOV)))
    end
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
    thisX=C{1}/200+abs(minX)+1;
    thisY=C{2}/200+abs(minY)+1;
    zombiecoo=floor(globalzombiespots(globalzombiespots(:,3)==(thisFOV),1:2));%transform the image to adjust camera orientation within mosaic
    zombiecoo(:,1)=2049-zombiecoo(:,1);
    zombiecoo(:,2)=2049-zombiecoo(:,2);
    zombieindeces=sub2ind([2048 2048],zombiecoo(:,1),zombiecoo(:,2));
    rotatedindeces=reverserotatedtemplate(zombieindeces);
    [zombiecoo(:,1),zombiecoo(:,2)]=ind2sub([2048 2048],rotatedindeces);
    zombiecoo(:,1)=floor((thisY-1)*1843.2+1+zombiecoo(:,1));%convert coordinates of decoded guide (zombie) in a FOV into global mosaic coordinates
    zombiecoo(:,2)=floor((thisX-1)*1843.2+1+zombiecoo(:,2));
    globalzombiespots(globalzombiespots(:,3)==(thisFOV),end-2:end)=[zombiecoo,rotatedindeces];
end
globalzombiespots(globalzombiespots(:,end-2)>0,end)=sub2ind([size(mymosaic,1) size(mymosaic,2)],globalzombiespots(globalzombiespots(:,end-2)>0,end-2),globalzombiespots(globalzombiespots(:,end-2)>0,end-1));
writematrix(globalzombiespots,fullfile(matlabtempfiles,'globalzombiespotsround3filteredQComemore.csv'))
%%
load(fullfile(matlabtempfiles,'finalmask.mat'))
clear globalmerfishspots merlinmask
cellsstats=regionprops(finalmask,'Centroid','PixelIdxList','Eccentricity','Area','Circularity','Extent');
clear mymosaic
countableZombie=zeros(size(cellsstats,1),9+78);

disp(strcat(num2str(size(cellsstats)),' objects detected'))
% globalmerfishspots(:,end+1)=zeros(size(globalmerfishspots,1),1);
globalzombiespots(:,end+1)=zeros(size(globalzombiespots,1),1);
for i=1:size(cellsstats,1)
    mypixels=cellsstats(i).PixelIdxList;%check for each cell if there is a decoded guide in it
%     positions=ismember(globalmerfishspots(:,10),mypixels);
    positionszombie=ismember(globalzombiespots(:,9),mypixels);
    %sum(double(positions));
%     globalmerfishspots(positions,11)=i;
    globalzombiespots(positionszombie,10)=i;%add labels to each guide with the ID of the cell they are part of
end
writematrix(globalzombiespots,fullfile(matlabtempfiles,'globalzombiespotswithcellIDbacktomaxround3filteredQComemore.csv'))
% writematrix(globalmerfishspots,fullfile(matlabtempfiles,'globalmerfishspotswithcellIDbacktomax.csv'))


for thiscell=1:size(cellsstats,1)%build the count table
    if cellsstats(thiscell).Area>0 
        countableZombie(thiscell,1)=1;%thisFOV;%#index of that cell
        countableZombie(thiscell,2)=thiscell;%#index of that cell
        countableZombie(thiscell,3)=cellsstats(thiscell).Area;
        countableZombie(thiscell,4)=cellsstats(thiscell).Eccentricity;
        countableZombie(thiscell,5)=cellsstats(thiscell).Circularity;
        countableZombie(thiscell,6)=cellsstats(thiscell).PixelIdxList(1);
        countableZombie(thiscell,7)=cellsstats(thiscell).Extent(1);
        mytable2=globalzombiespots(globalzombiespots(:,10)==thiscell,:);
        if size(mytable2,1)>0
             countableZombie(thiscell,8)=size(mytable2,1);
             for j=1:size(mytable2,1)
                 countableZombie(thiscell,8+mytable2(j,4))=1;
             end
        end
    end
end
writematrix(countableZombie,fullfile(matlabtempfiles,'zombietableonmosaicbacktomaxround3filteredQComemore.csv'))
