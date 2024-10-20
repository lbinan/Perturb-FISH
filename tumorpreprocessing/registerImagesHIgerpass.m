function registeredImage=registerImagesHIgerpass(imageref,thisimage)
fixed=imageref(:,:,28);%probably better with 29
data=fixed;
histFcn = @(x) histogram(x, 'Normalization', 'pdf');
h=histFcn(data);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
   ccum(i)=sum(h.BinCounts(1:i));
end          
x=[];
thre=0.993;
idx=find(ccum>thre*ccum(end));  
bw=data>h.BinEdges(idx(1));     
bw2=imfill(bw,'holes');
bw3=bwareaopen(bw2,10); 
Centroid = regionprops3(bw3, 'Centroid');
cfixed=Centroid.Centroid;
bw4=bw3;
%%
moving=thisimage(:,:,22);%and this one probably should stay 22
data=moving;
h=histFcn(data);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
   ccum(i)=sum(h.BinCounts(1:i));
end          
x=[];
thre=0.995;
idx=find(ccum>thre*ccum(end));  
bw=data>h.BinEdges(idx(1));     
bw2=imfill(bw,'holes');
bw3=bwareaopen(bw2,10); 
Centroid = regionprops3(bw3, 'Centroid');
cmove=Centroid.Centroid;
moving=[];
if sum(sum(bw3))==0 | sum(sum(bw4))==0
    registeredImage=thisimage;
else
%% shifting the beads
disp('Start warping');
cmovet=cmove;
registeredImage=zeros(2048,2048,23);
problem=0;
for trial=1:10
    if size(cmovet,2)>1 & size(cfixed,2)>1
     dist=pdist2(cmovet(:,1:2),cfixed(:,1:2));
    [M,I]=min(dist);
    [temp,ids]=sort(M);
    dev=diff(temp);
    idd=(dev<0.3*median(dev));
    temp1=cfixed(ids(idd),:);
    temp2=cmovet(I(ids(idd)),:);
    co=temp1-temp2;
coordi=[median(co(:,1)),median(co(:,2)),median(co(:,3))];    
cmovet=[cmovet(:,1)+coordi(1,1),cmovet(:,2)+coordi(1,2),cmovet(:,3)+coordi(1,3)];
    else
        problem=1;
    end
end
shift(1)=mean(cmovet(:,1)-cmove(:,1));
shift(2)=mean(cmovet(:,2)-cmove(:,2));
shift(3)=mean(cmovet(:,3)-cmove(:,3));
shift(isinf(shift)|isnan(shift))=0;
shiftpre=shift;
%% redo for weird shifts:
if max(abs(shift))>200 | sum(sum(bw4))<2500 | sum(sum(bw3))<2500
    disp('redoing weird one')
fixed=imageref(:,:,28);
data=fixed;
histFcn = @(x) histogram(x, 'Normalization', 'pdf');
h=histFcn(data);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
   ccum(i)=sum(h.BinCounts(1:i));
end          
x=[];
thre=0.95;
idx=find(ccum>thre*ccum(end));  
bw=data>h.BinEdges(idx(1));     
bw2=imfill(bw,'holes');
bw3=bwareaopen(bw2,10); 
Centroid = regionprops3(bw3, 'Centroid');
cfixed=Centroid.Centroid;
bw4=bw3;
%%
moving=thisimage(:,:,22);
data=moving;
h=histFcn(data);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
   ccum(i)=sum(h.BinCounts(1:i));
end          
x=[];
thre=0.96;
idx=find(ccum>thre*ccum(end));  
bw=data>h.BinEdges(idx(1));     
bw2=imfill(bw,'holes');
bw3=bwareaopen(bw2,10); 
Centroid = regionprops3(bw3, 'Centroid');
cmove=Centroid.Centroid;
moving=[];
%% shifting the beads
disp('Start warping');
cmovet=cmove;
registeredImage=zeros(2048,2048,23);
for trial=1:10
if size(cmovet,2)>1 & size(cfixed,2)>1
     dist=pdist2(cmovet(:,1:2),cfixed(:,1:2));
    [M,I]=min(dist);
    [temp,ids]=sort(M);
    dev=diff(temp);
    idd=(dev<0.3*median(dev));
    temp1=cfixed(ids(idd),:);
    temp2=cmovet(I(ids(idd)),:);
    co=temp1-temp2;
coordi=[median(co(:,1)),median(co(:,2)),median(co(:,3))];    
cmovet=[cmovet(:,1)+coordi(1,1),cmovet(:,2)+coordi(1,2),cmovet(:,3)+coordi(1,3)];
end
end
shift(1)=mean(cmovet(:,1)-cmove(:,1));
shift(2)=mean(cmovet(:,2)-cmove(:,2));
shift(3)=mean(cmovet(:,3)-cmove(:,3));
shift(isinf(shift)|isnan(shift))=0;
shiftpre=shift;
disp('weird one became')
shift
end
if max(abs(shift))>200 | sum(sum(bw4))<2500 | sum(sum(bw3))<2500
    disp('redoing weird one AGAIN')
fixed=imageref(:,:,28);
data=fixed;
histFcn = @(x) histogram(x, 'Normalization', 'pdf');
h=histFcn(data);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
   ccum(i)=sum(h.BinCounts(1:i));
end          
x=[];
thre=0.95;
idx=find(ccum>thre*ccum(end));  
bw=data>h.BinEdges(idx(1));     
bw2=imfill(bw,'holes');
bw3=bwareaopen(bw2,10); 
Centroid = regionprops3(bw3, 'Centroid');
cfixed=Centroid.Centroid;
bw4=bw3;
%%
moving=thisimage(:,:,22);
data=moving;
h=histFcn(data);
ccum=zeros(1,h.NumBins);
for i=1:h.NumBins
   ccum(i)=sum(h.BinCounts(1:i));
end          
x=[];
thre=0.95;
idx=find(ccum>thre*ccum(end));  
bw=data>h.BinEdges(idx(1));     
bw2=imfill(bw,'holes');
bw3=bwareaopen(bw2,10); 
Centroid = regionprops3(bw3, 'Centroid');
cmove=Centroid.Centroid;
moving=[];
%% shifting the beads
disp('Start warping');
cmovet=cmove;
registeredImage=zeros(2048,2048,23);
for trial=1:10
if size(cmovet,2)>1 & size(cfixed,2)>1
     dist=pdist2(cmovet(:,1:2),cfixed(:,1:2));
    [M,I]=min(dist);
    [temp,ids]=sort(M);
    dev=diff(temp);
    idd=(dev<0.3*median(dev));
    temp1=cfixed(ids(idd),:);
    temp2=cmovet(I(ids(idd)),:);
    co=temp1-temp2;
coordi=[median(co(:,1)),median(co(:,2)),median(co(:,3))];    
cmovet=[cmovet(:,1)+coordi(1,1),cmovet(:,2)+coordi(1,2),cmovet(:,3)+coordi(1,3)];
end
end
shift(1)=mean(cmovet(:,1)-cmove(:,1));
shift(2)=mean(cmovet(:,2)-cmove(:,2));
shift(3)=mean(cmovet(:,3)-cmove(:,3));
shift(isinf(shift)|isnan(shift))=0;
shiftpre=shift;
disp('weird one became')
shift
end
%% finish

[shift(1),shift(2)]
if problem==0
for i=1:size(thisimage,3)
registeredImage(:,:,i)=imtranslate(thisimage(:,:,i),[shift(1),shift(2)]);
end
else 
    registeredImage(:,:,i)=zeros(size(fixed));
end
end