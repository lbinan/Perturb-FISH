density=sum(pooledguides);
%here we split the data in subsets of increasing numbers of cells
%effects and p values are computed the "naive way" here to get a feeling of
%the data. they need to be properly computed using FR-perturb (see python
%script), using the tables that this script saves.
validation=[];
test6=[];
test12=[];
test23=[];
test40=[];
test50=[];
zombievalidation=[];
zombie6=[];
zombie12=[];
zombie23=[];
zombie40=[];
zombie50=[];
effects6=[];
effects12=[];
effects23=[];
effects40=[];
effects50=[];
effectsvalidations=[];
pcontrol=[];
p6=[];
p12=[];
p23=[];
p40=[];
p50=[];
for i=1:35
    numberofcells=density(i);
    thesecells=finalmerfish(pooledguides(:,i)==1,9:end);
    Z=zeros(1,36);
    Z(i)=1;
    for j=floor(numberofcells/2):numberofcells
        validation=[validation;thesecells(j,:)];
        zombievalidation=[zombievalidation;Z];
    end
    effectsvalidations=[effectsvalidations;log(mean(thesecells)./meancontrols)];
    [h,p]=ttest2(thesecells,controlcells(:,9:end));
    pcontrol=[pcontrol;p];
    for j=1:6
        test6=[test6;thesecells(j,:)];
        zombie6=[zombie6;Z];
    end
    effects6=[effects6;log(mean(thesecells(1:6,:))./meancontrols)];
    for j=1:12
        test12=[test12;thesecells(j,:)];
        zombie12=[zombie12;Z];
    end
    effects12=[effects12;log(mean(thesecells(1:12,:))./meancontrols)];
    for j=1:23
        test23=[test23;thesecells(j,:)];
        zombie23=[zombie23;Z];
    end
     for j=1:min(density(i),40)
        test40=[test40;thesecells(j,:)];
        zombie40=[zombie40;Z];
     end
    for j=1:min(density(i),50)
        test50=[test50;thesecells(j,:)];
        zombie50=[zombie50;Z];
    end
    
    effects23=[effects23;log(mean(thesecells(1:23,:))./meancontrols)];
    [h,p]=ttest2(thesecells(1:23,:),controlcells(:,9:end));
    p23=[p23;p];
    [h,p]=ttest2(thesecells(1:12,:),controlcells(:,9:end));
    p12=[p12;p];
    [h,p]=ttest2(thesecells(1:6,:),controlcells(:,9:end));
    p6=[p6;p];
    if density(i)>39
        effects40=[effects40;log(mean(thesecells(1:40,:))./meancontrols)];
    else
         effects40=[effects40;zeros(1,130)];
    end
    if density(i)>49
        effects50=[effects50;log(mean(thesecells(1:50,:))./meancontrols)];
    else
        effects50=[effects50;zeros(1,130)];
    end
        [h,p]=ttest2(thesecells(1:min(50,size(thesecells,1)),:),controlcells(:,9:end));
    p50=[p50;p];
    [h,p]=ttest2(thesecells(1:min(40,size(thesecells,1)),:),controlcells(:,9:end));
    p40=[p40;p];
end
controlcells=finalmerfish(pooledguides(:,37)==1|pooledguides(:,38)==1|pooledguides(:,39)==1,:);
meancontrols=mean(controlcells(:,9:end));
numberofcontrols=size(controlcells,1);
controlZ=zeros(numberofcontrols,36);
controlZ(:,36)=1;
zombievalidation=[zombievalidation;controlZ];
zombie6=[zombie6;controlZ];
zombie12=[zombie12;controlZ];
zombie23=[zombie23;controlZ];
validation=[validation;controlcells(:,9:end)];
test6=[test6;controlcells(:,9:end)];
test12=[test12;controlcells(:,9:end)];
test23=[test23;controlcells(:,9:end)];
test40=[test40;controlcells(:,9:end)];
test50=[test50;controlcells(:,9:end)];
mypath='\\helium\broad_clearylab\Users\Loic\thp1homemadezombie_1\numberofcellstests'
writematrix(validation,fullfile(mypath,'validationagain.csv'))
writematrix(test6,fullfile(mypath,'test6again.csv'))
writematrix(test12,fullfile(mypath,'test12again.csv'))
writematrix(test23,fullfile(mypath,'test23again.csv'))
writematrix(test40,fullfile(mypath,'test40again.csv'))
writematrix(test50,fullfile(mypath,'test50again.csv'))
writematrix(zombievalidation',fullfile(mypath,'zombievalidation.csv'))
writematrix(zombie6',fullfile(mypath,'zombie6.csv'))
writematrix(zombie12',fullfile(mypath,'zombie12.csv'))
writematrix(zombie23',fullfile(mypath,'zombie23.csv'))
writematrix(zombie40',fullfile(mypath,'zombie40.csv'))
writematrix(zombie50',fullfile(mypath,'zombie50.csv'))
zombie40=[zombie40;controlZ];
zombie50=[zombie50;controlZ];
writematrix(zombie40',fullfile(mypath,'zombie40.csv'))
writematrix(zombie50',fullfile(mypath,'zombie50.csv'))


%%
effects12cells=4*effects12cellsLFCs(:,1:35)';
effects6cells=4*effects6cellsLFCs(:,1:35)';
effects23cells=4*effects23cellsLFCs(:,1:35)';
effects40cells=4*effects40cellsLFCs(:,1:35)';
effects50cells=4*effects50cellsLFCs(:,1:35)';
effectsvalidationcells=4*effectsvalidationcellsLFCs(:,1:35)';
validationsQ=effectsvalidationcellsqvals(:,1:35)';
Q23cellsQ=effects23cellsqvals(:,1:35)';
figure,
for i=1:35
    if sum(validationsQ(i,:)<0.05)>0
        scatter(effects6cells(i,validationsQ(i,:)<0.05),effectsvalidationcells(i,validationsQ(i,:)<0.05),'.')
        hold on
    end
end
figure,
for i=1:35
    if sum(validationsQ(i,:)<0.05)>0
        scatter(effects12cells(i,validationsQ(i,:)<0.05),effectsvalidationcells(i,validationsQ(i,:)<0.05),'.')
        hold on
    end
end
figure,
for i=1:35
    if sum(pcontrol(i,:)<0.05)>0%&sum(effects50(i,:))~=0
        scatter(effectsvalidations(i,pcontrol(i,:)<0.05),effects6(i,pcontrol(i,:)<0.05),'.')
        hold on
    end
end
title('naive results,6 cells')

figure,heatmap(genelist,singleguides,effects50,'Colormap',parula(18))
figure,heatmap(genelist,singleguides,effects40,'Colormap',parula(18))
figure,heatmap(genelist,singleguides,effects23,'Colormap',parula(18))

figure,
for i=1:35
    if (sum(validationsQ(i,:)<0.05)>0)&(sum(abs(effects50cells(i,:))>0.8)>0)%&sum(effects50(i,:))~=0
        scatter(effectsvalidationcells(i,(validationsQ(i,:)<0.05)&(abs(effects50cells(i,:))>0.8)),effects50cells(i,(validationsQ(i,:)<0.05)&(abs(effects50cells(i,:))>0.8)),'.')
        hold on
    end
end
title('FR results,40 cells')

figure,heatmap(genelist,singleguides,effects12,'Colormap',parula(18))
figure,heatmap(genelist,singleguides,effects6,'Colormap',parula(18))



figure,
for i=1:35
     if (sum(validationsQ(i,:)<0.05)>0)%&(sum(abs(effects50cells(i,:))>0.8)>0)%&sum(effects50(i,:))~=0
        scatter(effectsvalidationcells(i,(validationsQ(i,:)<0.05)),effects40(i,(validationsQ(i,:)<0.05)),'.')
        hold on
    end
end
title('40 cells FR agains 40 cells naive')


