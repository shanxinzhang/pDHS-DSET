function output_label=pDHS_DSET(inputfile,outputfile,species,pwd_to_python, pwd_to_bioanalysis_feature)
%This function is used to prediciton of DNase Hypersensitive Sites in plant
%genome by using DS evidence based theory.
%Output: 
%output_label: vector of the prediction results, 1(DHSs), -1(Non
%DHSs), we also write the results to the outputfile at the same time
%Input: 
%inputfile: String type variety, the name of the input file (Fasta format), containing the sequences to be  predicted
%outpurfile: String type variety, the name of the outpur file
%species: plant species, 1-A.thaliana 2-rice
%pwd_to_python: absolute dirction of  python 
%pwd_to_bioanalysis_feature: 
%absolute dirction of  feature.py file in Bioseq-analysis program, please download the Bioseq-Analysis (http://bioinformatics.hitsz.edu.cn/BioSeq-Analysis/download) at first, and install it.

%Useage:
%output_label=pDHS_DSET('test.fas','test_result.txt',1,pwd_to_python, pwd_to_bioanalysis)

%Author: Shanxin Zhang
%Insitute: Jiangnan University
%E-mail:shanxinzhang@jiangnan.edu.cn



%0. get sequence information
[Header, ~] = fastaread(inputfile);
%1 generate features   
%1.1 generate kmer features
disp('Step1. Generating Various Features!');
disp('     1.1 Generating Kmer Features!');
     for k=1:5
         command=[pwd_to_python, ' ', pwd_to_bioanalysis_feature,'  ', inputfile, ' DNA -method Kmer -k ', num2str(k),' -r 0 -f tab -out DHSs_kmer_',num2str(k),'.txt'];
         system(command);
     end
    kmer_data=[];
    for l=1:5
        kmer_file=['DHSs_kmer_',num2str(l),'.txt'];
        data=importdata(kmer_file);
        delete(kmer_file);
        kmer_data=[kmer_data,data];
    end
    labels=ones(size(kmer_data,1),1);
     clear k command l kmer_file  data;
    %1.2 generate revckmer features
    disp('     1.2 Generating RevcKmer Features!');
     for k=1:5
     command=[pwd_to_python, ' ', pwd_to_bioanalysis_feature,'  ', inputfile, ' DNA -method Kmer -k ', num2str(k),' -r 1 -f tab -out DHSs_reckmer_',num2str(k),'.txt'];
     system(command);
     end
    %
    reckmer_data=[];
    for l=1:5
        reckmer_file=['DHSs_reckmer_',num2str(l),'.txt'];
        data=importdata(reckmer_file);
        delete(reckmer_file);
        reckmer_data=[reckmer_data,data];
    end
% [~,~,reckmer_decision_values]=libsvmpredict(labels,reckmer_data,reckmer_model,'-b 1 -q');  
  clear k command l reckmer_file  data;
    %1.3 mismatch features
    disp('     1.3 Generating Mismatch Profile Features!');
     for k=3:5
         command=[pwd_to_python, ' ', pwd_to_bioanalysis_feature,'  ', inputfile, ' DNA -method Mismatch -k ', num2str(k),' -m 1 -f tab -out DHSs_mismatch_',num2str(k),'.txt'];
         system(command);
     end
        mismatch_data=[];
    for l=3:5
        mismatch_file=['DHSs_mismatch_',num2str(l),'.txt'];
        data=importdata(mismatch_file);
        delete(mismatch_file);
        mismatch_data=[mismatch_data,data./repmat(sum(data,2),1,size(data,2))];
    end
     clear k command l mismatch_file  data;
     disp('     1.4 Generating PseDNC Features!');
    l=4;
   w=0.1;
   command=[pwd_to_python, ' ', pwd_to_bioanalysis_feature,'  ', inputfile, ' DNA -method PseDNC -lamada ', num2str(l), ' -w ', num2str(w), ' -f tab -out DHSs_psednc.txt'];
    system(command);
    psednc_file='DHSs_psednc.txt';
    psednc_data=importdata(psednc_file);
    delete(psednc_file);
    clear k command l psednc_file ;

disp('Step 2. Start Predicting  Various Features.');
if species==1 % prediction for A.thana 
    load('./TAIR_SVM_model/kmer_model_file.mat');
    load('./TAIR_SVM_model/reckmer_model_file.mat');
    load('./TAIR_SVM_model/mismatch_model_file.mat');
    load('./TAIR_SVM_model/psednc_model_file.mat');
else %Prediction for rice
    load('./TIGR_SVM_model/kmer_model_file.mat');
    load('./TIGR_SVM_model/reckmer_model_file.mat');
    load('./TIGR_SVM_model/mismatch_model_file.mat');
    load('./TIGR_SVM_model/psednc_model_file.mat');
end
    %generate kmer features prediction result
    disp('     2.1 Predicting Kmer Model.');
    [~,~,kmer_decision_values]=libsvmpredict(labels,kmer_data,kmer_model,'-b 1 -q'); 
    %1.2 generate revckmer features prediction result
    disp('     2.2 Predicting RevcKmer Model.');
    [~,~,reckmer_decision_values]=libsvmpredict(labels,reckmer_data,reckmer_model,'-b 1 -q');  
    disp('     2.3 Predicting Mismatch Profile Model.');
    [~,~,mismatch_decision_values]=libsvmpredict(labels,mismatch_data,mismatch_model,'-b 1 -q');  
   disp('     2.4 Predicting PseDNC Model.');
   [~,~,psednc_decision_values]=libsvmpredict(labels,psednc_data,psednc_model,'-b 1 -q');  
   disp('Step 3. Start Combining Prediction Results Based on DSET.');
     InputScores=[kmer_decision_values(:,1),reckmer_decision_values(:,1),mismatch_decision_values(:,1),psednc_decision_values(:,1)];
     [~,output_label]=DSET(InputScores,0.9);
     disp('Step 4. Write Prediction Results to Outputfile.');
     fid=fopen(outputfile,'w');
     for i=1:length(output_label)
         if length(output_label)==1
             if output_label(i)==1
                fprintf(fid, 'The %s is DHSs sequence.\n', Header);
             else
                 fprintf(fid,'The %s is Non DHSs sequence.\n', Header);
             end
         else
             if output_label(i)==1
                 fprintf(fid,'The %s is DHSs sequence.\n', Header{i});
             else
                 fprintf(fid,'The %s is Non DHSs sequence.\n', Header{i});
             end
         end
     end
     fclose(fid);
     disp('Step 5. Finish.');
end