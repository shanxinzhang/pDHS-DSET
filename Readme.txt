The code is written in MATLAB in Windows Enviroment. The python and pse-in-one tool is also employed. please modify the absolute direction of these two tools in the DHS_DSET.m file at first.Then, 
to predict the DHSs, please type the following command in MATLAB environment:
issuccess=DHS_DSET(inputfile,outputfile),
pleat note that the inputfile is fasta format.
For example:
issuccess=DHS_predict(¡®test.fasta¡¯,¡¯test_out.txt¡¯)

