%here as indication, the script actually used in the manuscript is
%checkingsubsamplecorrelation
%all the redoXcellsLFCs and redoXcellsqvals are the outputs of FR-perturb
%with X cells. they need to all be loaded to run this script
figure,heatmap(genelist,singleguides,redo50cellsLFCs(:,1:35)','Colormap',parula(18))
title('50 cells, FR')
figure,heatmap(genelist,singleguides,redo40cellsLFCs(:,1:35)','Colormap',parula(18))
title('40 cells, FR')
figure,heatmap(genelist,singleguides,redo23cellsLFCs(:,1:35)','Colormap',parula(18))
title('23 cells, FR')
figure,heatmap(genelist,singleguides,redo12cellsLFCs(:,1:35)','Colormap',parula(18))
title('12 cells, FR')
figure,heatmap(genelist,singleguides,redo6cellsLFCs(:,1:35)','Colormap',parula(18))
title('6 cells, FR')
figure,heatmap(genelist,singleguides,redovalidationcellsLFCs(:,1:35)','Colormap',parula(18))
title('validation cells, FR')
figure,heatmap(genelist,singleguides,redo50cellsLFCs(:,1:35)','Colormap',parula(18))
title('50 cells, FR,testnorm')

redovalidation=redovalidationcellsLFCs(:,1:35)';
redo50=redo50cellsLFCs(:,1:35)';
redo40=redo40cellsLFCs(:,1:35)';
redo23=redo23cellsLFCs(:,1:35)';
redo12=redo12cellsLFCs(:,1:35)';
redo6=redo6cellsLFCs(:,1:35)';
qvals50=redo50cellsqvals(:,1:35)';
qvals40=redo40cellsqvals(:,1:35)';
qvals23=redo23cellsqvals(:,1:35)';
qvals12=redo12cellsqvals(:,1:35)';
qvals6=redo6cellsqvals(:,1:35)';
qvalsvalidations=redovalidationcellsqvals(:,1:35)';

for i=1:35
    scatter(redovalidation(i,qvalsvalidations(i,:)<0.05),redo50(i,qvalsvalidations(i,:)<0.05),'.')
    hold on
end

redo50norm=testnormon50cellsLFCs(:,1:35)';
figure,heatmap(genelist,singleguides,redo50norm,'Colormap',parula(18))
title('50 cells, FR,testnorm')
redo50rank=testnormonpercent50cellsLFCs(:,1:35)';
figure,heatmap(genelist,singleguides,redo50rank,'Colormap',parula(18))
title('50 cells, FR,test norm0.1, rank10')


%% building precision and recal
flattenedQvalsValidation=reshape(qvalsvalidations,[(size(qvalsvalidations,1)*size(qvalsvalidations,2)) 1]);
flattenedQvals50cells=reshape(qvals50,[(size(qvals50,1)*size(qvals50,2)) 1]);
flattenedeffectsValidation=reshape(redovalidation,[(size(qvalsvalidations,1)*size(qvalsvalidations,2)) 1]);
flattenedeffects50cells=reshape(redo50,[(size(qvals50,1)*size(qvals50,2)) 1]);

valuestosort=[flattenedQvalsValidation,flattenedQvals50cells,flattenedeffectsValidation,flattenedeffects50cells];
sortedvalues=sortrows(valuestosort);
recall=[];
for i=1:size(sortedvalues,1)
%     if sortedvalues(i,1)<0.05
        recall=[recall;sum(sortedvalues(1:i,2)<0.05)/i];
%     end
end

%% building precision and recal
flattenedpvalsValidation=reshape(pcontrol,[(size(pcontrol,1)*size(pcontrol,2)) 1]);
flattenedpvals50cells=reshape(p50,[(size(p50,1)*size(p50,2)) 1]);
flattenedeffectsValidation=reshape(effectsvalidations,[(size(pcontrol,1)*size(pcontrol,2)) 1]);
flattenedeffects50cells=reshape(effects50,[(size(qvals50,1)*size(qvals50,2)) 1]);

valuestosort=[flattenedpvals50cells,flattenedpvalsValidation,flattenedeffectsValidation,flattenedeffects50cells];
sortedvalues=sortrows(valuestosort);
recall=[];
for i=1:size(sortedvalues,1)
        croppedtruth=sortedvalues(1:i,2)<0.05;
        croppedtest=sortedvalues(1:i,1)<0.05;
        TP=sum(croppedtruth>0&croppedtest>0);
        FP=sum(croppedtruth==0&croppedtest>0);
        FN=sum(croppedtruth>0&croppedtest==0);
        recall(i)=TP/(TP+FN);
        precision(i)=TP/(TP+FP);
end
scatter(recall,precision,'.')

flattenedpvalsValidation=reshape(pcontrol,[(size(pcontrol,1)*size(pcontrol,2)) 1]);
flattenedpvals40cells=reshape(p40,[(size(p50,1)*size(p50,2)) 1]);
flattenedeffectsValidation=reshape(effectsvalidations,[(size(pcontrol,1)*size(pcontrol,2)) 1]);
flattenedeffects40cells=reshape(effects40,[(size(qvals50,1)*size(qvals50,2)) 1]);

valuestosort=[flattenedpvals40cells,flattenedpvalsValidation,flattenedeffectsValidation,flattenedeffects40cells];
sortedvalues=sortrows(valuestosort);
recall=[];
for i=1:size(sortedvalues,1)
    croppedtruth=sortedvalues(i:end,2)<0.05;
    croppedtest=sortedvalues(i:end,1)<0.05;
    TP=sum(croppedtruth>0&croppedtest>0);
    FP=sum(croppedtruth==0&croppedtest>0);
    FN=sum(croppedtruth>0&croppedtest==0);
    recall(i)=TP/(TP+FN);
    precision(i)=TP/(TP+FP);
end
scatter(precision,recall,'.')

