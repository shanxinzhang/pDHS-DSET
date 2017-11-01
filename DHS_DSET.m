function issuccess=DHS_DSET(inputfile,outputfile)

% Prediction DHSs sequences from inputfile.
%   DHS_DSET(inputfile,outputfile) predict sequences from inputfile, and write results to outputfile
%  issucess=DHS_DSET(inputfile,outputfile), predict sequences from inputfile, and write results to outputfile, if success, return 1, else return 0
%   DHS_DSET(inputfile) predict sequences from inputfile, and write results to 'pred_result.txt' file.
%   Notes:
%   The inputfile is a fasta format file.
%   Examples:
% issucess=DHS_DSET('test.fasta','test_out.txt')
% Shanxin Zhang
% shanxinzhang@jiangnan.edu.cn


if nargin<1
    issuccess=0;
    disp('No input files found');
    return;
elseif nargin<2
    outputfile='pred_result.txt';
end

%please modify the direction of the python and pse-in-one tool in your own
%computer
 command=['C:\Python27\python.exe C:\Python27\Pse-in-One-1.0.3\Pse-in-One\acc.py -e user_indices.txt -f tab -l +1 -lag 1 ',inputfile,' DHSs_dac_1.txt DNA DAC'];
 system(command);
 command=['C:\Python27\python.exe C:\Python27\Pse-in-One-1.0.3\Pse-in-One\kmer.py -f tab -l +1 -r 1 -k 2 ', inputfile, ' DHSs_reckmer_2.txt DNA'];
 system(command);
 
 [head,~]=fastaread(inputfile);
 load ('model.mat');
 revckmer_data=importdata('DHSs_reckmer_2.txt');
 dac_data=importdata('DHSs_dac_1.txt');
 test_label=ones(size(revckmer_data,1),1);
 [~,~,revckmer_deci]=libsvmpredict(test_label,revckmer_data,revckmer_model,' -b 1 -q');  
  test_label=ones(size(dac_data,1),1);
 [~,~,dac_deci]=libsvmpredict(test_label,dac_data,dac_model,' -b 1 -q');  
  [~,pred_label]=DSET([revckmer_deci(:,1),dac_deci(:,1)]);
         fid=fopen(outputfile,'w');
        for i=1:length(head)
            if pred_label(i)==1
                fprintf(fid,'%s is predicted as DHSs.\n',head{i});
            else
                fprintf(fid,'%s is predicted as Non DHSs.\n',head{i});
            end
        end
        fclose(fid);
        clear all;
        issuccess=1;
end

