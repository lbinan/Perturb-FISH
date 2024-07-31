% neighbors=neighbors(filteredagainzombie(:,2),filteredagainzombie(:,2));

cellswith6neighbors=filteredagainmerfish(sum(neighbors,2)==6,9:end);
cellswith5neighbors=filteredagainmerfish(sum(neighbors,2)==5,9:end);
cellswith4neighbors=filteredagainmerfish(sum(neighbors,2)==4,9:end);
cellswith3neighbors=filteredagainmerfish(sum(neighbors,2)==3,9:end);
cellswith2neighbors=filteredagainmerfish(sum(neighbors,2)==2,9:end);
cellswith1neighbors=filteredagainmerfish(sum(neighbors,2)==1,9:end);
cellswith0neighbors=filteredagainmerfish(sum(neighbors,2)==0,9:end);


controlcellswith6neighbors=mean(filteredagainmerfish(sum(neighbors,2)==6,9:end),1);
controlcellswith5neighbors=mean(filteredagainmerfish(sum(neighbors,2)==5,9:end),1);
controlcellswith4neighbors=mean(filteredagainmerfish(sum(neighbors,2)==4,9:end),1);
controlcellswith3neighbors=mean(filteredagainmerfish(sum(neighbors,2)==3,9:end),1);
controlcellswith2neighbors=mean(filteredagainmerfish(sum(neighbors,2)==2,9:end),1);
controlcellswith1neighbors=mean(filteredagainmerfish(sum(neighbors,2)==1,9:end),1);
controlcellswith0neighbors=mean(filteredagainmerfish(sum(neighbors,2)==0,9:end),1);



bar([5718,6959,5608,3713,2064,1013])
title ("exact number of neighbors")
ylabel("cell count")

plottable=[log2(controlcellswith1neighbors./controlcellswith0neighbors);log2(controlcellswith2neighbors./controlcellswith0neighbors);log2(controlcellswith3neighbors./controlcellswith0neighbors);log2(controlcellswith4neighbors./controlcellswith0neighbors);log2(controlcellswith5neighbors./controlcellswith0neighbors);log2(controlcellswith6neighbors./controlcellswith0neighbors)];
figure, heatmap(genelist,["one neighbor","two neighbors","three neighbors", "four neighbors","five neighbors","six neighbors"],plottable,'Colormap',parula(18))
title('log2 effect of neighbors, regardless of guides, control have No neighbor')

[h,p1]=ttest2(cellswith1neighbors,cellswith0neighbors);
[h,p2]=ttest2(cellswith2neighbors,cellswith0neighbors);
[h,p3]=ttest2(cellswith3neighbors,cellswith0neighbors);
[h,p4]=ttest2(cellswith4neighbors,cellswith0neighbors);
[h,p5]=ttest2(cellswith5neighbors,cellswith0neighbors);
[h,p6]=ttest2(cellswith6neighbors,cellswith0neighbors);
plottable=[p1;p2;p3;p4;p5;p6];
figure, heatmap(genelist,["one neighbor","two neighbors","three neighbors", "four neighbors","five neighbors","six neighbors"],-log2(plottable),'Colormap',parula(18))
title('-log2 p values of neighbors, regardless of guides, control have No neighbor, all "significant"')

plottable=plottable<-log2(0.05);
figure, heatmap(genelist,["one neighbor","two neighbors","three neighbors", "four neighbors","five neighbors","six neighbors"],double(plottable),'Colormap',parula(18))

unperturbedcellsIDx=filteredagainzombie(:,9:82);
unperturbedcellsIDx=sum(unperturbedcellsIDx(:,71:74),2)>0&sum(unperturbedcellsIDx(:,1:70),2)<1;
sum(unperturbedcellsIDx)

controlcellswith6neighbors=filteredagainmerfish(unperturbedcellsIDx&sum(neighbors,2)==6,9:end);
controlcellswith5neighbors=filteredagainmerfish(unperturbedcellsIDx&sum(neighbors,2)==5,9:end);
controlcellswith4neighbors=filteredagainmerfish(unperturbedcellsIDx&sum(neighbors,2)==4,9:end);
controlcellswith3neighbors=filteredagainmerfish(unperturbedcellsIDx&sum(neighbors,2)==3,9:end);
controlcellswith2neighbors=filteredagainmerfish(unperturbedcellsIDx&sum(neighbors,2)==2,9:end);
controlcellswith1neighbors=filteredagainmerfish(unperturbedcellsIDx&sum(neighbors,2)==1,9:end);
controlcellswith0neighbors=filteredagainmerfish(unperturbedcellsIDx&sum(neighbors,2)==0,9:end);

guideID=[1:74];
guidesinallcells=guideID.*filteredagainzombie(:,9:82);
guidesinallcells=round(guidesinallcells/2);
guidesinallcells(guidesinallcells>35)=36;
%%
neighborsofControlcellswith2neighbors=neighbors(unperturbedcellsIDx&sum(neighbors,2)==2,:);
guidesinnighborsofcellswith2=zeros(size(neighborsofControlcellswith2neighbors,1),2);
for i=1:size(neighborsofControlcellswith2neighbors,1)
    theseneighbors=find(neighborsofControlcellswith2neighbors(i,:));
    if size(find(guidesinallcells(theseneighbors(1),:)),2)==1
        guidesinnighborsofcellswith2(i,1)=guidesinallcells(theseneighbors(1),find(guidesinallcells(theseneighbors(1),:)));
    else
        guidesinnighborsofcellswith2(i,1)=37;
    end
    if size(find(guidesinallcells(theseneighbors(2),:)),2)==1
        guidesinnighborsofcellswith2(i,2)=guidesinallcells(theseneighbors(2),find(guidesinallcells(theseneighbors(2),:)));
    else
        guidesinnighborsofcellswith2(i,2)=37;
    end
end
%%
neighborsofControlcellswith3neighbors=neighbors(unperturbedcellsIDx&sum(neighbors,2)==3,:);
guidesinnighborsofcellswith3=zeros(size(neighborsofControlcellswith3neighbors,1),3);
for i=1:size(neighborsofControlcellswith3neighbors,1)
    theseneighbors=find(neighborsofControlcellswith3neighbors(i,:));
    if size(find(guidesinallcells(theseneighbors(1),:)),2)==1
        guidesinnighborsofcellswith3(i,1)=guidesinallcells(theseneighbors(1),find(guidesinallcells(theseneighbors(1),:)));
    else
        guidesinnighborsofcellswith3(i,1)=37;
    end
    if size(find(guidesinallcells(theseneighbors(2),:)),2)==1
        guidesinnighborsofcellswith3(i,2)=guidesinallcells(theseneighbors(2),find(guidesinallcells(theseneighbors(2),:)));
    else
        guidesinnighborsofcellswith3(i,2)=37;
    end
    if size(find(guidesinallcells(theseneighbors(3),:)),2)==1
        guidesinnighborsofcellswith3(i,3)=guidesinallcells(theseneighbors(3),find(guidesinallcells(theseneighbors(3),:)));
    else
        guidesinnighborsofcellswith3(i,3)=37;
    end
end

%%
neighborsofControlcellswith4neighbors=neighbors(unperturbedcellsIDx&sum(neighbors,2)==4,:);
guidesinnighborsofcellswith4=zeros(size(neighborsofControlcellswith4neighbors,1),4);
for i=1:size(neighborsofControlcellswith4neighbors,1)
    theseneighbors=find(neighborsofControlcellswith4neighbors(i,:));
    if size(find(guidesinallcells(theseneighbors(1),:)),2)==1
        guidesinnighborsofcellswith4(i,1)=guidesinallcells(theseneighbors(1),find(guidesinallcells(theseneighbors(1),:)));
    else
        guidesinnighborsofcellswith4(i,1)=37;
    end
    if size(find(guidesinallcells(theseneighbors(2),:)),2)==1
    guidesinnighborsofcellswith4(i,2)=guidesinallcells(theseneighbors(2),find(guidesinallcells(theseneighbors(2),:)));
     else
        guidesinnighborsofcellswith4(i,2)=37;
    end
     if size(find(guidesinallcells(theseneighbors(3),:)),2)==1
    guidesinnighborsofcellswith4(i,3)=guidesinallcells(theseneighbors(3),find(guidesinallcells(theseneighbors(3),:)));
     else
        guidesinnighborsofcellswith4(i,3)=37;
     end
         if size(find(guidesinallcells(theseneighbors(4),:)),2)==1
    guidesinnighborsofcellswith4(i,4)=guidesinallcells(theseneighbors(4),find(guidesinallcells(theseneighbors(4),:)));
     else
        guidesinnighborsofcellswith4(i,4)=37;
    end
end
%%
neighborsofControlcellswith1neighbors=neighbors(unperturbedcellsIDx&sum(neighbors,2)==1,:);
guidesinnighborsofcellswith1=zeros(size(neighborsofControlcellswith4neighbors,1),1);
for i=1:size(neighborsofControlcellswith1neighbors,1)
    theseneighbors=find(neighborsofControlcellswith1neighbors(i,:));
    if size(find(guidesinallcells(theseneighbors(1),:)),2)==1
        guidesinnighborsofcellswith1(i,1)=guidesinallcells(theseneighbors(1),find(guidesinallcells(theseneighbors(1),:)));
    else
        guidesinnighborsofcellswith1(i,1)=37;
    end
end
%%
tableplot=[];
for i=1:37
    tableplot(i)=sum(guidesinnighborsofcellswith1==i);
end
bar(tableplot)
title("of all non perturbed cells that have one neighbor only, distribution of guides in that neighbor")
bar(tableplot(1:36))
title("of all non perturbed cells that have one neighbor only, distribution of guides in that neighbor")

tableplot=[];
for i=1:37
    tableplot(i)=sum(guidesinnighborsofcellswith2(:,1)==i|guidesinnighborsofcellswith2(:,2)==i);
end
bar(tableplot)
title("of all non perturbed cells that have TWO neighbors, distribution of cells having One neighbor with guide...")
bar(tableplot(1:36))
title("of all non perturbed cells that have TWO neighbors, distribution of cells having One neighbor with guide...")

tableplot=[];
for i=1:37
    tableplot(i)=sum(guidesinnighborsofcellswith2(:,1)==i&guidesinnighborsofcellswith2(:,2)==i);
end
bar(tableplot)
title("of all non perturbed cells that have TWO neighbors, distribution of cells having both neighbors with guide...")
bar(tableplot(1:36))
title("of all non perturbed cells that have TWO neighbors, distribution of cells having both neighbors with guide...")


tableplot=[];
for i=1:37
    tableplot(i)=sum(guidesinnighborsofcellswith3(:,1)==i|guidesinnighborsofcellswith3(:,2)==i|guidesinnighborsofcellswith3(:,2)==i);
end
bar(tableplot)
title("of all non perturbed cells that have THREE neighbors, distribution of cells having One neighbor with guide...")
bar(tableplot(1:36))
title("of all non perturbed cells that have THREE neighbors, distribution of cells having One neighbor with guide...")

tableplot=[];
for i=1:37
    tableplot(i)=sum(guidesinnighborsofcellswith3(:,1)==i&guidesinnighborsofcellswith3(:,2)==i&guidesinnighborsofcellswith3(:,3)==i);
end
bar(tableplot)
title("of all non perturbed cells that have THREE neighbors, distribution of cells having THREE neighbors with guide...")

%% effects from one guide in cells with 2 neighbors
controlcellswith2neighbors
keepeffects=ones(35,1);
effects=zeros(35,130);
for i=1:35
    testcellsIDx=guidesinnighborsofcellswith2(:,1)==i|guidesinnighborsofcellswith2(:,2)==i;
    if sum(testcellsIDx)<9
        keepeffects(i)=0;
    else
    testcounts=controlcellswith2neighbors(testcellsIDx,:);
    effects(i,:)=log2(mean(testcounts,1)./mean(controlcellswith2neighbors,1));
    [h,p]=ttest2(testcounts,controlcellswith2neighbors);
    pvaluestable(i,:)=p;
    end
end
figure, heatmap(genelist,singleguides(logical(keepeffects)),effects(logical(keepeffects),:),'Colormap',parula(18))
title('looking ONLY at non perturbed cells with 2 neighbors, effect of having a neighbor with guide X, showing only showing tests that occur 9 or more times')
figure, heatmap(genelist,singleguides(logical(keepeffects)),-log2(pvaluestable(logical(keepeffects),:)),'Colormap',parula(18))
title('looking ONLY at non perturbed cells with 2 neighbors, -log2(pvalues) of having a neighbor with guide X, showing only showing tests that occur 11 or more times')

table2cluster=effects(logical(keepeffects),:);
table2cluster(isinf(table2cluster))=0;
keptguides=singleguides(logical(keepeffects));
for i=1:size(keptguides,1)
    guidelabels{i}=keptguides(i);
end
cgo = clustergram(table2cluster,'Standardize','Row','Colormap',redbluecmap)
% set(cgo,'RowLabels',cellstr(keptguides),'ColumnLabels',cellstr(genelist),'Linkage','complete','Dendrogram',5)
set(cgo,'RowLabels',cellstr(keptguides),'ColumnLabels',cellstr(genelist))
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 10)
%%
normalizedtonumberofneighbors=cellswith0neighbors./mean(controlcellswith0neighbors,1);
normalizedtonumberofneighbors=[];
normalizedtonumberofneighbors=[normalizedtonumberofneighbors;controlcellswith1neighbors./mean(controlcellswith1neighbors,1);controlcellswith2neighbors./mean(controlcellswith2neighbors,1)];
normalizedtonumberofneighbors=[normalizedtonumberofneighbors;controlcellswith3neighbors./mean(controlcellswith3neighbors,1);controlcellswith4neighbors./mean(controlcellswith4neighbors,1)];
designfornormalizedtoneighbornumber=[];
designfornormalizedtoneighbornumber=[guidesinnighborsofcellswith1,zeros(size(guidesinnighborsofcellswith1,1),3)];
designfornormalizedtoneighbornumber=[designfornormalizedtoneighbornumber;[guidesinnighborsofcellswith2,zeros(size(guidesinnighborsofcellswith2,1),2)]];
designfornormalizedtoneighbornumber=[designfornormalizedtoneighbornumber;[guidesinnighborsofcellswith3,zeros(size(guidesinnighborsofcellswith3,1),1)]];
designfornormalizedtoneighbornumber=[designfornormalizedtoneighbornumber;[guidesinnighborsofcellswith4]];


keepeffects=ones(35,1);
effects=zeros(35,130);
for i=1:35
    testcellsIDx=designfornormalizedtoneighbornumber(:,1)==i|designfornormalizedtoneighbornumber(:,2)==i|designfornormalizedtoneighbornumber(:,3)==i|designfornormalizedtoneighbornumber(:,4)==i;
    if sum(testcellsIDx)<20
        keepeffects(i)=0;
    else
    testcounts=normalizedtonumberofneighbors(testcellsIDx,:);
    effects(i,:)=log2(mean(testcounts,1));
    [h,p]=ttest2(testcounts,controlcellswith4neighbors);
    pvaluestable(i,:)=p;
    end
end
figure, heatmap(genelist,singleguides(logical(keepeffects)),effects(logical(keepeffects),:),'Colormap',parula(18))
title('looking at effects of aving a neighbor with each of these guides, effects (appearing at least 20 times) compared to average in cells with same total number of neighbors, mean on all cells, across possible number of neighbors')
figure, heatmap(genelist([52,54]),singleguides([6,21,34]),effects([6,21,34],[52,54]),'Colormap',parula(18))
title('looking at effects of aving a neighbor with each of these guides, effects (appearing at least 20 times) compared to average in cells with same total number of neighbors, mean on all cells, across possible number of neighbors')
table2plot=effects(logical(keepeffects),:);
table2plot(isinf(table2plot))=0;
cgo = clustergram(table2plot,'Standardize','Row','Colormap',redbluecmap)
set(cgo,'RowLabels',cellstr(singleguides(logical(keepeffects))),'ColumnLabels',cellstr(genelist),'Linkage','complete','Dendrogram',5)
set(cgo,'RowLabels',cellstr(singleguides(logical(keepeffects))),'ColumnLabels',cellstr(genelist))
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 11)
