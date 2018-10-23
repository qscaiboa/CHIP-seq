# CHIP-seq


curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
python get-pip.py
pip -V

pip install MACS2
pip install MACS=1.4.2


#download miniconda python3 64bit
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install multiqc
conda install fastqc
conda install bowtie
conda install samblaster
conda install samtools
conda install deeptools
conda install R
conda install -c bioconda -c conda-forge snakemake
conda install sratoolkit


git clone https://github.com/kundajelab/phantompeakqualtools
R

install.packages("snow", repos="http://cran.us.r-project.org")
install.packages("snowfall", repos="http://cran.us.r-project.org")
install.packages("bitops", repos="http://cran.us.r-project.org")
install.packages("caTools", repos="http://cran.us.r-project.org")
source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools",suppressUpdates=TRUE)
install.packages("./spp_1.14.tar.gz")

ROSE


ascp installation


/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518123/SRR2518123.sra ./

/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518124/SRR2518124.sra ./

/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518125/SRR2518125.sra ./

/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518126/SRR2518126.sra ./


(wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - pi.dk/3) | 


wget -O - pi.dk/3
bash
find *sra| parallel -j 4  fastq-dump {}
