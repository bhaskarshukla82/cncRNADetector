REQUIREMENTS FOR CNCRNADETECTOR
1.	Blast download https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/
2.	Vienna web services https://www.tbi.univie.ac.at/RNA/
3.	Cygwin https://cygwin.com/
4.	Windows sub system(wsl) - Ubuntu – Microsoft store - https://www.microsoft.com/en-in/p/ubuntu-on-windows/9nblggh4msv6#activetab=pivot:overviewtab
5.	LGC - https://ngdc.cncb.ac.cn/biocode/tools/4/releases/4
6.	PLEK - https://sourceforge.net/projects/plek/files/
7.	CPC2 - http://cpc2.gao-lab.org/download.php
8.	pataa https://ftp.ncbi.nlm.nih.gov/blast/db/ 
9.	swissprot https://ftp.ncbi.nlm.nih.gov/blast/db/ 
10.	mature https://www.mirbase.org/ftp/CURRENT/
11.	tRNA http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/
12.	rRNA http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/
STEPS FOR INSTALLATION
1.	Download latest version of blast or as used in cncrnadetector blast 2.12.0. 
2.	For users working on windows platform download ncbi-blast-2.12.0+-win64.exe  whereas for Linux download ncbi-blast-2.12.0+-x64-linux.tar.gz (Links provided above)
3.	Next step would involve the installation of blast. In windows by clicking on the downloaded application setup appears. Follow the instructions and blast will be downloaded.
4.	For Linux users the tar file created need to be extracted in the specified directory.
5.	Create a folder named db in the bin of blast.
6.	Next step would involve downloading of the databases from NCBI (Links provided above).
7.	Download pataa (patented protein sequences) database as pataa.tar.gz file.
8.	Extract the file in the db folder.
9.	Download swissprot database as swissprot.tar.gz
10.	Extract the swissprot file.
11.	Create a folder named swissprot inside db folder
12.	Open terminal in the bin of blast and type the following command
12.1.1.	makeblastdb -in <swissprot_file_location> -dbtype prot -out db\swissprot\swissprot
13.	Download mature database as mature.fa.gz 
14.	Extract the mature file.
15.	Create a folder named mature inside db folder
16.	Open terminal in the bin of blast and type the following command
16.1.1.	makeblastdb -in <mature_file_location> -dbtype nucl -out db\swissprot\swissprot
17.	Download rfam above from the files 1 and 5 (Links provided above). 
18.	Download Vienna web services (Links provided above).
19.	Run setup to install Vienna web services.
20.	Download LGC-1.0.tar.gz (Links provided above).
21.	Download PLEK.1.2.tar.gz (Links provided above).
22.	Download CPC2 (Links provided above).
23.	For windows user install Cygwin or use ubuntu as subsystem.
