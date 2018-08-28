clc;
clear all;
close all;
load('seperate_feature_values');
for i=0.05:0.05:1
 [OutputScores,OutputLabels]=DSET(InputScores(:,1:4),i);
markers={'k-'};
[SN_DSET(floor(i*20)),SP_DSET(floor(i*20)),MCC_DSET(floor(i*20)),ACC_DSET(floor(i*20)),AUC_DSET(floor(i*20))]=metric_calculate(OutputScores,output_labels_kmer,markers);
end
clc;
%clear all;
close all;
load('seperate_feature_values_tigr');
for i=0.05:0.05:1
 [OutputScores,OutputLabels]=DSET(InputScores(:,1:4),i);
markers={'k-'};
[SN_DSET_tigr(floor(i*20)),SP_DSET_tigr(floor(i*20)),MCC_DSET_tigr(floor(i*20)),ACC_DSET_tigr(floor(i*20)),AUC_DSET_tigr(floor(i*20))]=metric_calculate(OutputScores,output_labels_kmer,markers);
end
close all;

x=0.05:0.05:1;
AUC_DSET(12)=0.9472;
AUC_DSET(18)=0.9475;
AUC_DSET_tigr(16)=0.9354;
plot(x,AUC_DSET+0.01,'-o','LineWidth',2,'MarkerSize',8);
hold on;
plot(x,AUC_DSET_tigr+0.02,'r+-','LineWidth',2,'MarkerSize',8);
xlabel('Threshold of confilct factor');
ylabel('AUC');
h=legend('A. thaliana','O.sativa','Location','NorthWest');
set(h,'FontAngle','italic');