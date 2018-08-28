The code is written in MATLAB in Windows Enviroment. The python and Bioseq-Analysis tool is also employed. please install Bioseq-Analysis and input the absolute direction of these two tools in the pDHS_DSET.m file at first.Then, 
to predict the DHSs, please type the following command in MATLAB environment:
output_labels=DHS_DSET(inputfile,outputfile,species, pwd_to_python, pwd_to_bioanalysis),
pleat note that the inputfile is fasta format.
For example:
output_label=pDHS_DSET('test.fas','test_result.txt',1,'C:\Python27\python.exe', 'D:\BioSeq-Analysis\feature.py');

