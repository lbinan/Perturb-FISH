%below are the rows on which to find each of the chosen example genes in
%our counttables
IRAK1=4
TRAF6=33
RELA=24
IRF3=6
NFKB1=20
MAP3K7=17

IRF3=6
TLR4=31
CD14=1
%goodalldataLFCs is the output from FR-perturb run on merged data from
%samples 1 and 2, goodalldataqvals is the matching Q-values, also from
%FR-perturb
effects1=goodalldataLFCs(:,IRF3);
Q1=goodalldataqvals(:,IRF3);
effects2=goodalldatasample2LFCs(:,IRF3);
Q2=goodalldatasample2qvals(:,IRF3);
scatter(effects1(Q1<0.1),effects2(Q1<0.1),'.')

subplot(2,2,1)
effects1=goodalldataLFCs(:,MAP3K7);
Q1=goodalldataqvals(:,MAP3K7);
effects2=goodalldatasample2LFCs(:,MAP3K7);
Q2=goodalldatasample2qvals(:,MAP3K7);
scatter(effects1(Q1<0.1),effects2(Q1<0.1),'.')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from MAP3K7 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'])
xlabel('Sample 1')
ylabel('Sample2')
subplot(2,2,2)
effects1=goodalldataLFCs(:,IRAK1);
Q1=goodalldataqvals(:,IRAK1);
effects2=goodalldatasample2LFCs(:,IRAK1);
Q2=goodalldatasample2qvals(:,IRAK1);
scatter(effects1(Q1<0.1),effects2(Q1<0.1),'.')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from IRAK1 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'])
xlabel('Sample 1')
ylabel('Sample2')
subplot(2,2,3)
effects1=goodalldataLFCs(:,TRAF6);
Q1=goodalldataqvals(:,TRAF6);
effects2=goodalldatasample2LFCs(:,TRAF6);
Q2=goodalldatasample2qvals(:,TRAF6);
scatter(effects1(Q1<0.1),effects2(Q1<0.1),'.')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from TRAF6 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'])
xlabel('Sample 1')
ylabel('Sample2')
subplot(2,2,4)

scatter(flateffect1(flatQ1<0.1&flatQ2<0.1),flateffect2(flatQ1<0.1&flatQ2<0.1),'.')
a=corrcoef(flateffect1(flatQ1<0.1&flatQ2<0.1),flateffect2(flatQ1<0.1&flatQ2<0.1))
title(['All significant effects, correlation:',' ', num2str(round(100*a(1,2))),'%'])
xlabel('Sample 1')
ylabel('Sample2')
set(gcf,'color','w');

set(gcf,'color','w');


effects1=goodalldataLFCs(:,NFKB1);
Q1=goodalldataqvals(:,NFKB1);
effects2=goodalldatasample2LFCs(:,NFKB1);
Q2=goodalldatasample2qvals(:,NFKB1);
scatter(effects1(Q1<0.1),effects2(Q1<0.1),'.')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from NFKB1 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'])
xlabel('Sample 1')
ylabel('Sample2')
set(gcf,'color','w');

flateffect1=reshape(goodalldataLFCs(:,1:34),[size(goodalldataLFCs(:,1:34),1)*size(goodalldataLFCs(:,1:34),2),1]);
flateffect2=reshape(goodalldatasample2LFCs(:,1:34),[size(goodalldataLFCs(:,1:34),1)*size(goodalldataLFCs(:,1:34),2),1]);
flatQ1=reshape(goodalldataqvals(:,1:34),[size(goodalldataLFCs(:,1:34),1)*size(goodalldataLFCs(:,1:34),2),1]);
flatQ2=reshape(goodalldatasample2qvals(:,1:34),[size(goodalldataLFCs(:,1:34),1)*size(goodalldataLFCs(:,1:34),2),1]);
flatDoug=reshape(neweffectsfromDoug,[size(neweffectsfromDoug,1)*size(neweffectsfromDoug,2),1]);
flatDougQ=reshape(dougsQvalues,[size(dougsQvalues,1)*size(dougsQvalues,2),1]);
corrcoef(flateffect1(flatQ1<0.1&flatQ2<0.1),flateffect2(flatQ1<0.1&flatQ2<0.1))
flateffect1=reshape(goodalldataLFCspool(:,1:34),[size(goodalldataLFCs(:,1:34),1)*size(goodalldataLFCs(:,1:34),2),1]);
flatQ1=reshape(goodalldataqvalspool(:,1:34),[size(goodalldataLFCs(:,1:34),1)*size(goodalldataLFCs(:,1:34),2),1]);

scatter(flateffect1(flatQ1<0.1&flatQ2<0.1),flateffect2(flatQ1<0.1&flatQ2<0.1),'.')
a=corrcoef(flateffect1(flatQ1<0.1&flatQ2<0.1),flateffect2(flatQ1<0.1&flatQ2<0.1))
title(['All significant effects, correlation:',' ', num2str(round(100*a(1,2))),'%'])
xlabel('Sample 1')
ylabel('Sample2')
set(gcf,'color','w');


corrcoef(flateffect1(flatQ1<0.1&flatDougQ<0.1),flatDoug(flatQ1<0.1&flatDougQ<0.1))

scatter(flatDoug(flatQ1<0.1&flatDougQ<0.1),flateffect1(flatQ1<0.1&flatDougQ<0.1),25,'filled')
a=corrcoef(flateffect1(flatQ1<0.1&flatDougQ<0.1),flatDoug(flatQ1<0.1&flatDougQ<0.1))
title(['All significant effects, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 10)
xlabel('Perturb-seq', 'FontSize', 10)
ylabel('Perturb-FISH', 'FontSize', 10)
set(gcf,'color','w');


scatter(flatDoug(flatQ1<0.1&flatDougQ<0.1),flateffect1(flatQ1<0.1&flatDougQ<0.1),25,'filled')
a=corrcoef(flateffect1(flatQ1<0.1&flatDougQ<0.1),flatDoug(flatQ1<0.1&flatDougQ<0.1))
title(['All significant effects, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 10)
xlabel('Perturb-seq', 'FontSize', 10)
ylabel('Perturb-FISH', 'FontSize', 10)
set(gcf,'color','w');


effects1=goodalldataLFCs(:,IRF3);
Q1=goodalldataqvals(:,IRF3);
effects2=goodalldatasample2LFCs(:,IRF3);
Q2=goodalldatasample2qvals(:,IRF3);
scatter(effects1(Q1<0.1),effects2(Q1<0.1),'.')

subplot(2,2,1)
effects1=goodalldataLFCs(:,MAP3K7);
Q1=goodalldataqvals(:,MAP3K7);
effects2=goodalldatasample2LFCs(:,MAP3K7);
Q2=goodalldatasample2qvals(:,MAP3K7);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects2(Q1<0.1),effects1(Q1<0.1))
title(['Effects from MAP3K7 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Sample 1', 'FontSize', 12)
ylabel('Sample 2', 'FontSize', 12)
subplot(2,2,2)
effects1=goodalldataLFCs(:,IRAK1);
Q1=goodalldataqvals(:,IRAK1);
effects2=goodalldatasample2LFCs(:,IRAK1);
Q2=goodalldatasample2qvals(:,IRAK1);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from IRAK1 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Sample 1', 'FontSize', 12)
ylabel('Sample 2', 'FontSize', 12)
xlim([-0.6 0.6])
ylim([-0.3 0.5])

subplot(2,2,3)
effects1=goodalldataLFCs(:,TRAF6);
Q1=goodalldataqvals(:,TRAF6);
effects2=goodalldatasample2LFCs(:,TRAF6);
Q2=goodalldatasample2qvals(:,TRAF6);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from TRAF6 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Sample 1', 'FontSize', 12)
ylabel('Sample 2', 'FontSize', 12)
xlim([-0.6 0.6])
ylim([-0.5 0.6])

subplot(2,2,4)

scatter(flateffect2(flatQ1<0.1&flatQ2<0.1),flateffect1(flatQ1<0.1&flatQ2<0.1),25,'filled')
a=corrcoef(flateffect1(flatQ1<0.1&flatQ2<0.1),flateffect2(flatQ1<0.1&flatQ2<0.1))
title(['All significant effects, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Sample 1', 'FontSize', 12)
ylabel('Sample 2', 'FontSize', 12)
set(gcf,'color','w');



%%
subplot(2,3,1)
effects1=goodalldataLFCs(:,IRF3);
Q1=goodalldataqvals(:,IRF3);
effects2=goodalldatasample2LFCs(:,IRF3);
Q2=goodalldatasample2qvals(:,IRF3);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from IRF3 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 10)
xlabel('Sample 1', 'FontSize', 10)
ylabel('Sample 2', 'FontSize', 10)
xlim([-0.6 0.6])
ylim([-0.3 0.5])

subplot(2,3,2)
effects1=goodalldataLFCs(:,TLR4);
Q1=goodalldataqvals(:,TLR4);
effects2=goodalldatasample2LFCs(:,TLR4);
Q2=goodalldatasample2qvals(:,TLR4);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from TLR4 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 10)
xlabel('Sample 1', 'FontSize', 10)
ylabel('Sample 2', 'FontSize', 10)
xlim([-0.6 0.6])
ylim([-0.5 0.6])

subplot(2,3,3)

effects1=goodalldataLFCs(:,CD14);
Q1=goodalldataqvals(:,CD14);
effects2=goodalldatasample2LFCs(:,CD14);
Q2=goodalldatasample2qvals(:,CD14);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from CD14 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 10)
xlabel('Sample 1', 'FontSize', 10)
ylabel('Sample 2', 'FontSize', 10)
xlim([-0.6 0.6])
ylim([-0.5 0.6])
subplot(2,3,4)

effects1=goodalldataLFCspool(:,IRF3);
Q1=goodalldataqvalspool(:,IRF3);
effects2=neweffectsfromDoug(:,IRF3);
Q2=dougsQvalues(:,IRF3);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from IRF3 KO'], 'FontSize', 10)
xlabel('Perturb-seq', 'FontSize', 10)
ylabel('Perturb FISH', 'FontSize', 10)
xlim([-0.6 0.6])
ylim([-0.3 0.5])

subplot(2,3,5)
effects1=goodalldataLFCspool(:,TLR4);
Q1=goodalldataqvalspool(:,TLR4);
effects2=neweffectsfromDoug(:,TLR4);
Q2=dougsQvalues(:,TLR4);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from TLR4 KO'], 'FontSize', 10)
xlabel('Perturb-seq', 'FontSize', 10)
ylabel('Perturb FISH', 'FontSize', 10)
xlim([-0.6 0.6])
ylim([-0.5 0.6])

subplot(2,3,6)

effects1=goodalldataLFCspool(:,CD14);
Q1=goodalldataqvalspool(:,CD14);
effects2=neweffectsfromDoug(:,CD14);
Q2=dougsQvalues(:,CD14);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from CD14 KO'], 'FontSize', 10)
xlabel('Perturb-seq', 'FontSize', 10)
ylabel('Perturb FISH', 'FontSize', 10)
xlim([-0.6 0.6])
ylim([-0.5 0.6])

set(gcf,'color','w');

%%

subplot(2,2,1)
effects1=goodalldataLFCspool(:,MAP3K7);
Q1=goodalldataqvalspool(:,MAP3K7);
effects2=neweffectsfromDoug(:,MAP3K7);
Q2=dougsQvalues(:,MAP3K7);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects2(Q1<0.1),effects1(Q1<0.1))
title(['Effects from MAP3K7 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Perturb-seq', 'FontSize', 12)
ylabel('Perturb-FISH', 'FontSize', 12)
subplot(2,2,2)
effects1=goodalldataLFCspool(:,IRAK1);
Q1=goodalldataqvalspool(:,IRAK1);
effects2=neweffectsfromDoug(:,IRAK1);
Q2=dougsQvalues(:,IRAK1);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from IRAK1 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Perturb-seq', 'FontSize', 12)
ylabel('Perturb-FISH', 'FontSize', 12)
xlim([-0.6 0.6])
ylim([-0.3 0.5])

subplot(2,2,3)
effects1=goodalldataLFCspool(:,TRAF6);
Q1=goodalldataqvalspool(:,TRAF6);
effects2=neweffectsfromDoug(:,TRAF6);
Q2=dougsQvalues(:,TRAF6);
scatter(effects2(Q1<0.1),effects1(Q1<0.1),25,'filled')
a=corrcoef(effects1(Q1<0.1),effects2(Q1<0.1))
title(['Effects from TRAF6 KO, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Perturb-seq', 'FontSize', 12)
ylabel('Perturb-FISH', 'FontSize', 12)
xlim([-0.6 0.6])
ylim([-0.5 0.6])

subplot(2,2,4)

scatter(flateffect1(flatQ1<0.1&flatDougQ<0.1),flatDoug(flatQ1<0.1&flatDougQ<0.1),25,'filled')
a=corrcoef(flateffect1(flatQ1<0.1&flatDougQ<0.1),flateffect2(flatQ1<0.1&flatDougQ<0.1))
title(['All significant effects, correlation:',' ', num2str(round(100*a(1,2))),'%'], 'FontSize', 12)
xlabel('Perturb-seq', 'FontSize', 12)
ylabel('Perturb-FISH', 'FontSize', 12)
set(gcf,'color','w');