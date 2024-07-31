mypath='\\helium\broad_clearylab\Users\Loic\thp1homemadezombie_1\downsampling10_30'
allsubsets=["5","10","20","25","30","40","50","validation"];%these are all the downsampled datasets
validationeffects=readmatrix(fullfile(mypath, 'effectsfromvalidationcells_LFCs.csv'));%these are the ouputs of FR-perturb on the validation dataset
validationQ=readmatrix(fullfile(mypath, 'effectsfromvalidationcells_qvals.csv'));
validationeffects=validationeffects([1:30,32:end],2:end-1);validationQ=validationQ([1:30,32:end],2:end-1);%removing genes that are missing from our matched single cell data
V=reshape(validationeffects,[size(validationeffects,1)*size(validationeffects,2), 1]);
Q=reshape(validationQ,[size(validationeffects,1)*size(validationeffects,2), 1]);
correlations=[];figure,
for i=1:7
    subsets=allsubsets(i);
    test=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_LFCs.csv')));%load FR outpus for each subsampled dataset
    test=test([1:30,32:end],2:end-1);
    testQ=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_qvals.csv')));
    testQ=testQ([1:30,32:end],2:end-1);
    testeffect=reshape(test,[size(validationeffects,1)*size(validationeffects,2), 1]);
    testQ=reshape(testQ,[size(validationeffects,1)*size(validationeffects,2), 1]);
    a=corrcoef(test(Q<0.1&testQ<0.1),V(Q<0.1&testQ<0.1));a(isnan(a))=0;
    b=corrcoef(test(Q<0.1&testQ<1),V(Q<0.1&testQ<1));
    correlations=[correlations;[i,b(1,2),a(1,2),sum(Q<0.1&testQ<0.1),sum(Q<0.1)]];
    subplot(2,4,i)
    if sum(Q<0.1&testQ<0.1)>0
        scatter(V(Q<0.1&testQ<1),testeffect(Q<0.1&testQ<1),60,'.')%plot those significant in validation
    end
    hold on
    if sum(Q<0.1&testQ<0.1)>0
        scatter(V(Q<0.1&testQ<0.1),testeffect(Q<0.1&testQ<0.1),80,'.r')%plot those significant in both
    end   
    ylabel(['Effects from', ' ',convertStringsToChars(subsets),' cells'])
    xlabel('Effects from validation set')
end
sgtitle('Scatter plots of effect sizes (LFCs) from different data subsets') 
set(gcf,'color','w');

plot([5,10,20,25,30,40,50],correlations(:,2),'-*','LineWidth',2)
hold on
plot([5,10,20,25,30,40,50],correlations(:,3),'-*','LineWidth',2)
plot([5,10,20,25,30,40,50],correlations(:,3),'-*','LineWidth',2)
plot([5,10,20,25,30,40,50],correlations(:,3),'-*','LineWidth',2)
plot([5,10,20,25,30,40,50],correlations(:,2),'-*','LineWidth',2)
set(gcf,'color','w');
ylim([0.4 1])
xlabel('Number of cells', 'FontSize', 14)
ylabel('Pearson correlation', 'FontSize', 14)
title('Correlation of effects with increasing numbers of cells per perturbation', 'FontSize', 14)



correlations=[];figure,
for i=1:8
    subsets=allsubsets(i);
    test=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_LFCs.csv')));
    test=test([1:30,32:end],2:end-1);
    testQ=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_qvals.csv')));
    testQ=testQ([1:30,32:end],2:end-1);
    testeffect=reshape(test,[size(validationeffects,1)*size(validationeffects,2), 1]);
    testQ=reshape(testQ,[size(validationeffects,1)*size(validationeffects,2), 1]);
    a=corrcoef(test(Q<0.1&testQ<0.1),V(Q<0.1&testQ<0.1));a(isnan(a))=0;
    b=corrcoef(test(Q<0.1&testQ<1),V(Q<0.1&testQ<1));
    correlations=[correlations;[i,b(1,2),a(1,2)]];
    subplot(2,3,i)
    if sum(Q<0.1&testQ<0.1)>0
        scatter(V(Q<0.1&testQ<1),testeffect(Q<0.1&testQ<1),'.')
    end
    title(strcat('Effects from', subsets,' cells'))
end
sgtitle('Effects significant in validation data subset') 
set(gcf,'color','w');
correlations=[["","correlation significant in validation","correlation significant in both"];correlations];
correlations=[["";allsubsets(1:6)'],correlations];
writematrix(correlations,fullfile(mypath, 'correlations.csv'))

allcellseffects=readmatrix(fullfile(mypath, 'goodalldata_LFCs.csv'));
allcellsQ=readmatrix(fullfile(mypath, 'goodalldata_qvals.csv'));
allcellseffects=allcellseffects(:,2:end-1);allcellsQ=allcellsQ(:,2:end-1);
allcellsVeffects=reshape(allcellseffects,[size(allcellseffects,1)*size(allcellseffects,2), 1]);
allcellsQ=reshape(allcellsQ,[size(allcellseffects,1)*size(allcellseffects,2), 1]);

validationeffects=readmatrix(fullfile(mypath, 'effectsfromvalidationcells_LFCs.csv'));
validationQ=readmatrix(fullfile(mypath, 'effectsfromvalidationcells_qvals.csv'));
validationeffects=validationeffects([1:30,32:88,90:end],[2:10,12:end-1]);
validationQ=validationQ([1:30,32:88,90:end],[2:10,12:end-1]);
validationQ=reshape(validationQ,[size(validationeffects,1)*size(validationeffects,2), 1]);
validationQ=reshape(validationQ,[size(validationeffects,1)*size(validationeffects,2), 1]);

truthsign=validationQ<0.1;
q=[0.001:0.001:1]
for i=1:7%build precision and recall tables
    subsets=allsubsets(i);
    test=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_LFCs.csv')));
    test=test([1:30,32:88,90:end],[2:10,12:end-1]);
    testQ=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_qvals.csv')));
    testQ=testQ([1:30,32:88,90:end],[2:10,12:end-1]);
    testeffect=reshape(test,[size(validationeffects,1)*size(validationeffects,2), 1]);
    testQ=reshape(testQ,[size(validationeffects,1)*size(validationeffects,2), 1]);
    R=[]
    P=[]
    for qID=1:size(q,2)
        P=[P;sum(truthsign&testQ<q(1+size(q,2)-qID))/sum(testQ<q(1+size(q,2)-qID))];
        R=[R;sum(truthsign&testQ<q(1+size(q,2)-qID))/sum(truthsign)];
    end
    line(P,R)
    hold on
end
% R=[]
% P=[]
% for qID=1:size(q,2)
%     P=[P;sum(truthsign&validationQ<q(1+size(q,2)-qID))/sum(validationQ<q(1+size(q,2)-qID))];
%     R=[R;sum(truthsign&validationQ<q(1+size(q,2)-qID))/sum(truthsign)];
% end
% 
% line(P,R)
% hold on
% R=[]
% P=[]
% for qID=1:size(q,2)
%     P=[P;sum(truthsign&allcellsQ<q(1+size(q,2)-qID))/sum(allcellsQ<q(1+size(q,2)-qID))];
%     R=[R;sum(truthsign&allcellsQ<q(1+size(q,2)-qID))/sum(truthsign)];
% end

line(P,R)
title('Precision recall')
xlabel('Precision')
ylabel('Recall')

%%
figure,
ax = axes();
hold(ax,'on');
truthsign=validationQ<0.1;
q=[0.001:0.001:1]
AUC=[];
for i=1:7%compute area under the curve
    subsets=allsubsets(i);
    test=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_LFCs.csv')));
    test=test([1:30,32:88,90:end],[2:10,12:end-1]);
    testQ=readmatrix(fullfile(mypath, strcat('effectsfrom',subsets,'cells_qvals.csv')));
    testQ=testQ([1:30,32:88,90:end],[2:10,12:end-1]);
    testeffect=reshape(test,[size(validationeffects,1)*size(validationeffects,2), 1]);
    testQ=reshape(testQ,[size(validationeffects,1)*size(validationeffects,2), 1]);
    R=[1]
    P=[0]
    for qID=1:size(q,2)
        P=[P;sum(truthsign&testQ<q(1+size(q,2)-qID))/sum(testQ<q(1+size(q,2)-qID))];
        R=[R;sum(truthsign&testQ<q(1+size(q,2)-qID))/sum(truthsign)];
    end
    R=[R;0]
    P=[P;1]
    ID=~isnan(R)&~isnan(P);
    R=R(ID);
    P=P(ID);

    % line(P,R)
    plot(P,R,'DisplayName',[strcat(subsets, " Cells")],'LineWidth',2)
    AUC=[AUC;trapz(P,R)];
    % hold on
end

% R=[]
% P=[]
% for qID=1:size(q,2)
%     P=[P;sum(truthsign&validationQ<q(1+size(q,2)-qID))/sum(validationQ<q(1+size(q,2)-qID))];
%     R=[R;sum(truthsign&validationQ<q(1+size(q,2)-qID))/sum(truthsign)];
% end
% 
% plot(P,R,'DisplayName',["half data"],'LineWidth',2)
% R=[]
% P=[]
% for qID=1:size(q,2)
%     P=[P;sum(truthsign&allcellsQ<q(1+size(q,2)-qID))/sum(allcellsQ<q(1+size(q,2)-qID))];
%     R=[R;sum(truthsign&allcellsQ<q(1+size(q,2)-qID))/sum(truthsign)];
% end

% plot(P,R,'DisplayName',["all data"])
title('Precision/recall curves for different numbers of cells', 'FontSize', 14)
xlabel('Precision', 'FontSize', 14)
ylabel('Recall', 'FontSize', 14)
hold(ax,'off');
legend();
set(gcf,'color','w');
legend('FontSize', 10);

figure,
plot([5,10,20,25,30,40,50],AUC,'-*','LineWidth',2,'DisplayName',['Perturb-FISH'])
hold on
plot([5,10,20,25,30,40,50],0.21*ones(1,7),'--','LineWidth',2,'DisplayName',['Perturb-seq approximate AUPRC for 50 cells'])

% title('Precision/recall curves for different numbers of cells', 'FontSize', 14)
xlabel('Number of cells', 'FontSize', 14)
ylabel('AUPRC', 'FontSize', 14)
ylim([0 0.4])
set(gcf,'color','w');
hold(ax,'off');
legend();
set(gcf,'color','w');
legend('FontSize', 10);



figure,
ax = axes();
hold(ax,'on');

plot([5,10,20,25,30,40,50],correlations(:,2),'-*','LineWidth',2,'DisplayName',['Perturb-FISH'])
hold on
plot([5,10,20,25,30,40,50],0.82*ones(1,7),'--','LineWidth',2,'DisplayName',['Perturb-seq approximate correlations for 50 cells'])
hold(ax,'off');
legend();

% title('Precision/recall curves for different numbers of cells', 'FontSize', 14)
xlabel('Number of cells', 'FontSize', 14)
ylabel('Pearson correlation', 'FontSize', 14)
ylim([0 1])
xlim([0 50])
set(gcf,'color','w');

set(gcf,'color','w');
legend('FontSize', 10);