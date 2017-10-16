#!/bin/bash

#Xubuntu install
sudo apt-get update
sudo apt-get upgrade
sudo apt-get dist-upgrade

#install build-essential (headers, etc.)
sudo apt-get install build-essential

#install cmake
sudo apt-get install cmake

#mount the extra HHD at startup
sudo mkdir /media/6tb_raid10
sudo cp /etc/fstab /etc/fstab.backup
echo "/dev/sdb1 /media/6tb_raid10 ext4 defaults 0 0 " | sudo tee -a /etc/fstab
sudo mount -a

#install system monitor
sudo apt install gnome-system-monitor

#install Disks
sudo apt install gnome-disk-utility

#install gedit
sudo apt install gedit

#install tree
sudo apt install tree

#install aptitude
sudo apt-get install aptitude

#update perl
sudo cpan "CPAN"
sudo cpan -r

#there qill be many prompts

#install java 8
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get update
sudo apt-get install oracle-java8-installer

#exfat file system support
sudo apt-get install exfat-utils exfat-fuse

#pigz
sudo apt-get install pigz

#Get latest R
echo "deb https://muug.ca/mirror/cran/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
sudo gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E084DAB9
sudo gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get update
sudo apt-get install r-base r-base-dev
#setup R local cran mirror
echo -e "local({
  r <- getOption(\"repos\")
  r[\"CRAN\"] <- \"https://muug.ca/mirror/cran/\"
  options(repos = r)
})
" > "${HOME}"/.Rprofile

#install RStudio
cd "${HOME}"/prog
sudo apt-get install libjpeg62-dev
sudo apt-get install libgstreamer-plugins-base1.0-dev
wget https://download1.rstudio.org/rstudio-1.0.143-amd64.deb
sudo dpkg -i *.deb
rm *.deb

#install GATK
cd "${HOME}"/prog
#must register to download -> https://www.broadinstitute.org/gatk/
mkdir gatk
cd gatk
tar xvjf GenomeAnalysisTK-3.6.tar.bz2
rm GenomeAnalysisTK-3.5.tar.bz2
chmod +x GenomeAnalysisTK.jar
echo -e "\n#GATK\nexport PATH=\$HOME/prog/gatk:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install htslib/samtools/bcftools
#install htslib
cd "${HOME}"/prog
wget https://github.com/samtools/htslib/releases/download/1.4/htslib-1.4.tar.bz2
tar xvjf htslib-1.4.tar.bz2
rm htslib-1.4.tar.bz2
cd htslib-1.4
./configure
make -j
sudo make install

#install samtools
cd "${HOME}"/prog
wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2
tar xvjf samtools-1.4.tar.bz2
rm samtools-1.4.tar.bz2
cd samtools-1.4
./configure
make -j
sudo make install

#install bcftools
cd "${HOME}"/prog
wget https://github.com/samtools/bcftools/releases/download/1.4/bcftools-1.4.tar.bz2
tar xvjf bcftools-1.4.tar.bz2
rm bcftools-1.4.tar.bz2
cd bcftools-1.4
./configure
make -j
sudo make install

#install sambamba
cd "${HOME}"/prog
wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2
tar xvjf sambamba_v0.6.6_linux.tar.bz2
rm sambamba_v0.6.6_linux.tar.bz2
mkdir sambamba-v0.6.6
mv sambamba_v0.6.6 sambamba-v0.6.6/sambamba
chmod +x sambamba-v0.6.6/sambamba
echo -e "\n#sambamba\nexport PATH=\$HOME/prog/sambamba-v0.6.6:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install bbmap
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/bbmap/BBMap_37.10.tar.gz
tar zxvf BBMap_37.10.tar.gz
rm BBMap_37.10.tar.gz
echo -e "\n#BBMap\nexport PATH=\$HOME/prog/bbmap:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install picard tools
cd "${HOME}"/prog
mkdir picard-tools
cd picard-tools
wget https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar
chmod +x picard-tools-2.5.0/*.jar
echo -e "\n#Picard-tools\nexport PATH=\$HOME/prog/picard-tools:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install bwa
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
tar xvjf bwa-0.7.15.tar.bz2
rm bwa-0.7.15.tar.bz2
cd bwa-0.7.15
make -j
sudo cp bwa /usr/local/bin

#install git
sudo apt-get install git
sudo apt-get install build-essential checkinstall
sudo apt-get install cvs subversion git-core mercurial

#install freebayes
cd "${HOME}"/prog
git clone --recursive https://github.com/ekg/freebayes.git
cd freebayes
make
sudo make install

#install Beagle
cd "${HOME}"/prog
mkdir beagle
cd beagle
wget https://faculty.washington.edu/browning/beagle/beagle.21Jan17.6cc.jar

#install RAxML
cd "${HOME}"/prog
git clone --recursive https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML
make -f Makefile.AVX2.PTHREADS.gcc
rm *.o
sudo cp raxmlHPC-PTHREADS-AVX2 /usr/local/bin

#install newick utilities
sudo apt-get update
sudo apt-get install libtool bison flex
cd "${HOME}"/prog
git clone --recursive https://github.com/tjunier/newick_utils.git
cd newick_utils
autoreconf -fi
./configure
make
sudo make install

#install rsvg-convert
sudo apt-get install librsvg2-bin

#install kSNP
sudo apt-get install tcsh
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/ksnp/kSNP3.021_Linux_package.zip
unzip kSNP3.021_Linux_package.zip
rm kSNP3.021_Linux_package.zip
rm -rf __MACOSX
mv kSNP3.021_Linux_package/kSNP3 .
#rm -rf kSNP3.021_Linux_package
cd kSNP3
sed -i 's%/usr/local/bin/kSNP3%/usr/local/bin%' kSNP3
sudo cp * /usr/local/bin
#get the source to modify is needed
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/ksnp/kSNP3.021_Source.zip
unzip kSNP3.021_Source.zip
rm kSNP3.021_Source.zip
rm -rf __MACOSX
#modify Kchooser from source and compile with "pp"
cd kSNP3.021_Source
sed -i 's/output_0/output/' Kchooser
#install pp
sudo apt install libpar-packer-perl
rm ../kSNP3/Kchooser
pp -o ../kSNP3/Kchooser Kchooser
sudo cp ../kSNP3/Kchooser /usr/local/bin


############## Prokka - Start ##############

#install blast+
cd "${HOME}"/prog
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz
rm ncbi-blast-2.6.0+-x64-linux.tar.gz
cd ncbi-blast-2.6.0+/bin
sudo cp * /usr/local/bin

#install expat
sudo apt-get install libexpat-dev

#install bioperl
sudo apt-get install cpanminus
sudo cpanm Bio::Perl

#install hmmer-2
cd "${HOME}"/prog
wget http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz
tar zxvf hmmer-2.3.2.tar.gz
rm hmmer-2.3.2.tar.gz
cd hmmer-2.3.2
./configure
make -j
sudo make install

#install hmmer-3
cd "${HOME}"/prog
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz
tar zxvf hmmer-3.1b2.tar.gz
rm hmmer-3.1b2.tar.gz
cd hmmer-3.1b2
./configure
make -j
sudo make install

#install aragorn
cd "${HOME}"/prog
wget http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/aragorn1.2.38.tgz
tar zxvf aragorn1.2.38.tgz
rm aragorn1.2.38.tgz
cd aragorn1.2.38
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.38.c
sudo cp aragorn /usr/local/bin

#install prodigal
cd "${HOME}"/prog
git clone --recursive https://github.com/hyattpd/Prodigal.git
sudo make install

#install tbl2asn
cd "${HOME}"/prog
mkdir tbl2asn-linux64
cd tbl2asn-linux64
wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
gunzip -d linux64.tbl2asn.gz
chmod +x linux64.tbl2asn
mv linux64.tbl2asn tbl2asn
sudo cp tbl2asn /usr/local/bin

#install GNU parallel
sudo apt-get install parallel

#install infernal
cd "${HOME}"/prog
wget http://eddylab.org/infernal/infernal-1.1.2.tar.gz
tar zxvf infernal-1.1.2.tar.gz
rm infernal-1.1.2.tar.gz
cd infernal-1.1.2
./configure
make -j
sudo make install

#install barrnap
cd "${HOME}"/prog
wget http://www.vicbioinformatics.com/barrnap-0.6.tar.gz
tar zxvf barrnap-0.6.tar.gz
rm barrnap-0.6.tar.gz
cd barrnap-0.6
echo -e "\n#barrnap\nexport PATH=\$HOME/prog/barrnap-0.6/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install MINCED
cd "${HOME}"/prog
git clone --recursive https://github.com/ctSkennerton/minced.git
cd minced
make -j
sudo cp minced /usr/local/bin
sudo cp minced.jar /usr/local/bin/

#install rnammer
# see http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer
cd "${HOME}"/prog
mkdir rnammer
mv rnammer-1.2.src.tar.Z rnammer/
cd rnammer
tar zxvf rnammer-1.2.src.tar.Z
rm rnammer-1.2.src.tar.Z
echo -e "\n#rnammer\nexport PATH=\$HOME/prog/rnammer:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
#set the vairiable "$INSTALL_PATH" to "/home/duceppemo/prog/rnammer"
#set the vairiable "$HMMSEARCH_BINARY" to "/home/duceppemo/prog/hmmer-2.3.2/src/hmmsearch"
sed -i "s%/usr/cbs/bio/src/rnammer-1.2%$PWD%" rnammer
sed -i "s%/usr/cbs/bio/bin/linux64/hmmsearch%$HOME/prog/hmmer-2.3.2/src/hmmsearch%" rnammer
# Remove the two instances of –cpu 1 (not the whole sentence, just "–cpu 1") in the lines containing "system sprintf('%s --cpu 1 --compat"
sed -i 's/ --cpu 1//' core-rnammer

#install signalP
#see http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp
cd "${HOME}"/prog
tar zxvf signalp-4.1f.Linux.tar.gz
rm signalp-4.1f.Linux.tar.gz
cd signalp-4.1
sed -i "s%/usr/cbs/bio/src/signalp-4.1%$HOME/prog/signalp-4.1%" signalp
echo -e "\n#Perl library\nexport PERL5LIB=\$HOME/prog/signalp-4.1/lib:$PERL5LIB" | tee -a "${HOME}"/.bashrc
echo -e "\n#signalP\nexport PATH=\$HOME/prog/signalp-4.1:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install some perl modules
#http://sourceforge.net/projects/swissknife/
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/swissknife/swissknife/1.73/Swissknife_1.73.tar.gz
tar zxvf Swissknife_1.73.tar.gz
rm Swissknife_1.73.tar.gz
cd Swissknife_1.73
perl Makefile.PL
make
make test
sudo make install

#install tmhmm
sudo apt-get install gnuplot
# see http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm
cd "${HOME}"/prog
tar zxvf tmhmm-2.0c.Linux.tar.gz
rm tmhmm-2.0c.Linux.tar.gz
cd tmhmm-2.0c/bin
sudo cp decodeanhmm.Linux_x86_64 /usr/local/bin/

#install Prokka
sudo apt-get install libdatetime-perl libxml-simple-perl
sudo cpan Bio::Root:Version
sudo cpan Parallel::Loops
sudo cpan Sys::CPU
cd "${HOME}"/prog
git clone --recursive https://github.com/tseemann/prokka.git
echo -e "\n#Prokka\nexport PATH=\$HOME/prog/prokka/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
ln -s /usr/local/bin/barrnap "${HOME}"/prog/prokka/binaries/linux/
#update kindom db
cd prokka/bin
prokka --setupdb
./prokka-build_kingdom_dbs
prokka --setupdb
#update HAMAP db
cd ../db/hmm
wget ftp://ftp.expasy.org/databases/hamap/rules_index.dat
wget ftp://ftp.expasy.org/databases/hamap/hamap_rules.dat
wget ftp://ftp.expasy.org/databases/hamap/hamap_seed_alignments.tar.gz
wget ftp://ftp.expasy.org/databases/hamap/hamap.prf.gz
tar zxvf hamap_seed_alignments.tar.gz
gunzip hamap.prf.gz
cp ../../bin/prokka-hamap_to_hmm ../../bin/prokka-hamap_to_hmm.backup
# make compatible with new file names
sed -i 's/alignment_id.dat/rules_index.dat/' ../../bin/prokka-hamap_to_hmm
sed -i 's/hamap_families.dat/hamap_rules.dat/' ../../bin/prokka-hamap_to_hmm
#make prokka-hamap_to_hmm parallel
sed -i 's/File::Temp qw(tempdir);/File::Temp qw(tempdir);\nuse Parallel::Loops;\nuse Sys::CPU;\n\nmy $cpus = Sys::CPU::cpu_count();\nmy $pl = Parallel::Loops->new($cpus);/' ../../bin/prokka-hamap_to_hmm
sed -i 's/for my $msa (@msa) {/$pl->foreach( \\@msa, sub {\n    my $msa = $_;/' ../../bin/prokka-hamap_to_hmm
sed -i '/product name/{n;n;s/\}/\}\)\;/;}' ../../bin/prokka-hamap_to_hmm
perl ../../bin/prokka-hamap_to_hmm --datadir=$PWD #or run my modified script that downloads the files
rm HAMAP.hmm.*
hmmpress HAMAP.hmm
prokka --setupdb

#Genus db
#Custom built bast protein db
cd ../genus
cp file.fasta ./
makeblastdb -in "${HOME}"/file.fasta \
    -dbtype "prot" \
    -input_type "fasta" \
    -title "B07-007_Chr" \
    -parse_seqids \
    -hash_index \
    -out ./Listeria


############## Prokka - End ##############

#install snp-sites
sudo apt-get install snp-sites

#install mothur
cd "${HOME}"/prog
wget https://github.com/mothur/mothur/releases/download/v1.39.5/Mothur.linux_64.zip
unzip Mothur.linux_64.zip
rm Mothur.linux_64.zip
rm -rf __MACOSX
echo -e "\n#Mothur\nexport PATH=\$HOME/prog/mothur:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install zlib
sudo apt-get install lib64z1-dev

#install readline
sudo apt-get install libreadline6-dev

#install bowtie2
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.1/bowtie2-2.3.1-linux-x86_64.zip
unzip bowtie2-2.3.1-linux-x86_64.zip
rm bowtie2-2.3.1-linux-x86_64.zip
cd bowtie2-2.3.1
sudo cp bowtie2* /usr/local/bin

#install Trinity
cd "${HOME}"/prog
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.4.0.tar.gz
tar zxvf Trinity-v2.4.0.tar.gz
rm Trinity-v2.4.0.tar.gz
cd trinityrnaseq-Trinity-v2.4.0
make -j
make plugins
#test install
cd sample_data/test_Trinity_Assembly/
./runMe.sh
echo -e "\n#Trinity\nexport PATH=\$HOME/prog/trinityrnaseq-2.2.0:\$HOME/prog/trinityrnaseq-2.2.0/utils:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install Transdecoder
cd "${HOME}"/prog
git clone --recursive https://github.com/TransDecoder/TransDecoder.git
cd TransDecoder
make -j
echo -e "\n#TransDecoder\nexport PATH=\$HOME/prog/TransDecoder:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install SQLite
cd "${HOME}"/prog
wget https://sqlite.org/2017/sqlite-autoconf-3180000.tar.gz
tar zxvf sqlite-autoconf-3180000.tar.gz
rm sqlite-autoconf-3180000.tar.gz
cd sqlite-autoconf-3180000
./configure
make -j
sudo make install

#install Trinotate
cd "${HOME}"/prog
wget https://github.com/Trinotate/Trinotate/archive/v3.0.2.tar.gz
tar zxvf v3.0.2.tar.gz
rm v3.0.2.tar.gz
cd Trinotate-3.0.2
echo -e "\n#Trinotate\nexport PATH=\$HOME/prog/Trinotate-3.0.2:\$PATH" | tee -a "${HOME}"/.bashrc
#Swissprot db
mkdir /media/6tb_raid10/db/swissprot
cd /media/6tb_raid10/db/swissprot
wget https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz
gunzip uniprot_sprot.pep.gz
makeblastdb -in uniprot_sprot.pep -dbtype prot
#Pfam-A db
cd /media/6tb_raid10/db/pfam-a
wget https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
#SQLite db
mkdir /media/6tb_raid10/db/trinotateSQLite
cd /media/6tb_raid10/db/trinotateSQLite
wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz" -O Trinotate.sqlite.gz
gunzip Trinotate.sqlite.gz

#install Kallisto
sudo apt-get install libhdf5-dev
cd "${HOME}"/prog
git clone --recursive https://github.com/pachterlab/kallisto.git
cd kallisto
cmake CMakeLists.txt 
make -j
sudo cp src/kallisto /usr/local/bin

#install sleuth
sudo apt-get install curl openssl
sudo R
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("pachterlab/sleuth")
library('sleuth')

#install RSEM
cd "${HOME}"/prog
git clone --recursive https://github.com/deweylab/RSEM.git
cd RSEM
make -j
sudo make install

#install DEseq2
sudo R
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")

#install fastQC
cd "${HOME}"/prog
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
rm fastqc_v0.11.5.zip
cd FastQC
chmod +x fastqc
echo -e "\n#FastQC\nexport PATH=\$HOME/prog/FastQC:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install QUAST
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/quast/quast-4.5.tar.gz
tar zxvf quast-4.5.tar.gz
rm quast-4.5.tar.gz
cd quast-4.5
sudo ./setup.py install_full
# sudo ./install_full.sh
# echo -e "\n#quast\nexport PATH=\$HOME/prog/quast-4.5/quast.py:\$PATH" | tee -a "${HOME}"/.bashrc
# source "${HOME}"/.bashrc

#install e-mem
# see here for source: http://www.csd.uwo.ca/~ilie/E-MEM/
cd "${HOME}"/prog
unzip e-mem.zip
rm -rf __MACOSX
make

#install LASER
sudo apt-get install libboost-dev
cd "${HOME}"/prog
wget http://www.csd.uwo.ca/~ilie/LASER/LASER.zip
unzip LASER.zip
rm -rf __MACOSX
cd LASER/libs/MUMmer3.23-linux
rm -rf e-mem-linux/
wget http://www.csd.uwo.ca/~ilie/E-MEM/e-mem.zip
mv e-mem_2 e-mem-linux/
cd ../..
python quast.py --test
python metaquast.py --test
######Wont compile.######

#install SWIG
cd "${HOME}"/prog
wget http://prdownloads.sourceforge.net/swig/swig-3.0.12.tar.gz
tar zxvf swig-3.0.12.tar.gz
rm swig-3.0.12.tar.gz
cd swig-3.0.12
./configure
make -j
sudo make install

#install Jellyfish
cd "${HOME}"/prog
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz
tar zxvf jellyfish-2.2.6.tar.gz
rm jellyfish-2.2.6.tar.gz
cd jellyfish-2.2.6
sudo apt-get install python-dev
./configure --enable-python-binding --enable-perl-binding
make -j
sudo make install

#install cd-hit
cd "${HOME}"/prog
git clone --recursive https://github.com/weizhongli/cdhit.git
cd cdhit
make -j
cd cd-hit-auxtools
make -j
echo -e "\n#cd-hit\nexport PATH=\$HOME/prog/cdhit:\$HOME/prog/cdhit/cd-hit-auxtools:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install Tassel/UNEAK
cd "${HOME}"/prog
git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git
cd tassel-5-standalone
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"
sed -i "s/-Xmx1536m/"$memJava"/" start_tassel.pl


#install Stacks
cd "${HOME}"/prog
wget http://catchenlab.life.illinois.edu/stacks/source/stacks-1.46.tar.gz
tar zxvf stacks-1.46.tar.gz
rm stacks-1.46.tar.gz
cd stacks-1.46
./configure
make -j
sudo make install

#install popoolation
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/popoolation/popoolation_1.2.2.zip
# wget https://downloads.sourceforge.net/project/popoolation2/popoolation2_1201.zip
unzip popoolation_1.2.2.zip
rm popoolation_1.2.2.zip
cd popoolation_1.2.2
for i in $(find $(pwd) -type f | grep ".pl");do chmod +x "$i";done
echo -e "\n#Perl library for popoolation\nexport PERL5LIB=\$HOME/prog/popoolation_1.2.2/Modules:$PERL5LIB" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install PLINK 1.9
cd "${HOME}"/prog
mkdir plink
cd plink
wget https://www.cog-genomics.org/static/bin/plink170330/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip
rm plink_linux_x86_64.zip
sudo cp plink /usr/local/bin
sudo cp prettify /usr/local/bin

#install PLINK 2
cd "${HOME}"/prog
mkdir plink2
cd plink2
wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20170420.zip
unzip plink2_linux_x86_64_20170420.zip
rm plink2_linux_x86_64_20170420.zip
sudo cp plink2 /usr/local/bin

#install jellyfish 1
cd "${HOME}"/prog
wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz
tar zxvf jellyfish-1.1.11.tar.gz
rm jellyfish-1.1.11.tar.gz
cd jellyfish-1.1.11
./configure --prefix=$PWD
make
make install
cd bin
mv jellyfish jellyfish1
sudo cp jellyfish1 /usr/local/bin

#install Kraken -> requires jellyfish v1
cd "${HOME}"/prog
wget https://ccb.jhu.edu/software/kraken/dl/kraken-0.10.5-beta.tgz
tar zxvf kraken-0.10.5-beta.tgz
rm kraken-0.10.5-beta.tgz
cd kraken-0.10.5-beta
./install_kraken.sh $PWD
echo -e "\n#kraken\nexport PATH=\$HOME/prog/kraken-0.10.5-beta:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
sed -i 's/jellyfish --version/jellyfish1 --version/' check_for_jellyfish.sh
sed -i 's/jellyfish count/jellyfish1 count/' build_kraken_db.sh
sed -i 's/jellyfish merge/jellyfish1 merge/' build_kraken_db.sh
cd srcipts
sed -i 's/jellyfish --version/jellyfish1 --version/' check_for_jellyfish.sh
sed -i 's/jellyfish count/jellyfish1 count/' build_kraken_db.sh
sed -i 's/jellyfish merge/jellyfish1 merge/' build_kraken_db.sh
#build custom db
#download taxonomy
mkdir -p /media/6tb_raid10/db/kraken/refseq_BV_old/taxonomy
cd /media/6tb_raid10/db/kraken/refseq_BV_old/taxonomy
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar zxvf taxdump.tar.gz
rm taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
pigz -d gi_taxid_nucl.dmp.gz
#download refseq all bacterial genomes (old) and viruses (latest)
mkdir -p /media/6tb_raid10/db/kraken/refseq_BV_old/refseq
cd media/6tb_raid10/db/kraken/refseq_BV_old/refseq
wget ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz
tar zxvf all.fna.tar.gz
rm all.fna.tar.gz
mkdir -p /media/6tb_raid10/db/kraken/refseq_BV_old/viruses_fna
cd media/6tb_raid10/db/kraken/refseq_BV_old/viruses_fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz #complete genome
tar zxvf all.fna.tar.gz
rm all.fna.tar.gz
mkdir -p /media/6tb_raid10/db/kraken/refseq_BV_old/viruses_fna
cd media/6tb_raid10/db/kraken/refseq_BV_old/viruses_fna 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.ffn.tar.gz #protein coding regions
tar zxvf all.ffn.tar.gz
rm all.ffn.tar.gz
#prepare library with fasta files
mkdir -p /media/6tb_raid10/db/kraken/refseq_BV_old/library
for i in $(find /media/6tb_raid10/db/kraken/refseq_BV_old -type f | grep "ffn\|fna"); do
    mv "$i" /media/6tb_raid10/db/kraken/refseq_BV_old/library
done
rm -rf /media/6tb_raid10/db/kraken/refseq_BV_old/refseq \
    /media/6tb_raid10/db/kraken/refseq_BV_old/viruses_ffn \
    /media/6tb_raid10/db/kraken/refseq_BV_old/viruses_fna
#build kraken db
kraken-build --build \
    --db /media/6tb_raid10/db/kraken/refseq_BV_old \
    --threads $(nproc) \
    --jellyfish-hash-size $(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100))
#clean unused files by kraken
kraken-build --clean \
    --db /media/6tb_raid10/db/kraken/refseq_BV_old \
    --threads $(nproc)


#install Krona
cd "${HOME}"/prog
wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
tar xvf KronaTools-2.7.tar
rm KronaTools-2.7.tar
cd KronaTools-2.7
mkdir /media/6tb_raid10/db/taxonomy
sudo ./install.pl --taxonomy /media/6tb_raid10/db/taxonomy
./updateTaxonomy.sh
./updateAccessions.sh

#install circos
#sudo apt-get install circos
sudo apt-get install libgd-dev

#circos -module
sudo cpan 'SVG'
sudo cpan 'YAML'
sudo cpan 'Config::General'
sudo cpan 'Font::TTF'
sudo cpan 'GD'
sudo cpan 'List::MoreUtils'
sudo cpan 'Test::Regexp'
sudo cpan 'Math::Bezier'
sudo cpan 'Math::Round'
sudo cpan 'Math::VecStat'
sudo cpan 'Params::Validate'
sudo cpan 'Readonly'
sudo cpan 'Regexp::Common'
sudo cpan 'Set::IntSpan'
sudo cpan 'Text::Format'
sudo cpan 'Statistics::Basic'

cd /bin
sudo ln -s /usr/bin/env env

cd "${HOME}"/prog
wget http://circos.ca/distribution/circos-0.69-4.tgz
tar zxvf circos-0.69-4.tgz
rm circos-0.69-4.tgz
cd circos-0.69-4
echo -e "\n#Circos\nexport PATH=\$HOME/prog/circos-0.69-4/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

 #install circos tools
cd "${HOME}"/prog/circos-0.69-4
wget http://circos.ca/distribution/circos-tools-0.22.tgz
tar zxvf circos-tools-0.22.tgz
rm circos-tools-0.22.tgz
cd circos-tools-0.22/tools/binlinks
./run

#install circos tutorials
cd "${HOME}"/prog/circos-0.69-4
wget http://circos.ca/distribution/circos-tutorials-0.67.tgz
tar zxvf circos-tutorials-0.67.tgz
rm circos-tutorials-0.67.tgz
#run generate tutorial image files
for i in $(find "$HOME/prog/circos-0.69-4/circos-tutorials-0.67/tutorials" -type d); do
    for j in $(find "$i" -type d);do
        [ -e "${j}"/circos.conf ] && circos -conf "${j}"/circos.conf
    done
done

#install mauve
cd "${HOME}"/prog
wget http://darlinglab.org/mauve/snapshots/2015/2015-02-13/linux-x64/mauve_linux_snapshot_2015-02-13.tar.gz
tar zxvf mauve_linux_snapshot_2015-02-13.tar.gz
rm mauve_linux_snapshot_2015-02-13.tar.gz
echo -e "\n#Mauve\nexport PATH=\$HOME/prog/mauve_snapshot_2015-02-13:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install smalt
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/smalt/smalt-0.7.6-static.tar.gz
tar zxvf smalt-0.7.6-static.tar.gz
rm smalt-0.7.6-static.tar.gz
cd smalt-0.7.6
./configure
make -j
sudo make install

#install REAPR
sudo cpan 'File::Basename'
sudo cpan 'File::Copy'
sudo cpan 'File::Spec'
sudo cpan 'File::Spec::Link'
sudo cpan 'Getopt::Long'
sudo cpan 'List::Util'
cd "${HOME}"/prog
wget ftp://ftp.sanger.ac.uk/pub/resources/software/reapr/Reapr_1.0.18.tar.gz
tar zxvf Reapr_1.0.18.tar.gz
rm Reapr_1.0.18.tar.gz
cd Reapr_1.0.18
./install.sh
# modify to use latest smalt and samtools instead of the ones provided
sed -i 's/samtools sort $raw_bam $raw_bam.sort/samtools sort -@ $options{n} -o $raw_bam.sort.bam $raw_bam/' task_smaltmap.pl
sed -i 's/samtools view -S -T/samtools view -@ $options{n} -S -T/' task_smaltmap.pl
sed -i 's/samtools view -H/samtools view -@ $options{n} -H/' task_smaltmap.pl
sed -i 's/$smalt index/smalt index/' task_smaltmap.pl
sed -i 's/$smalt sample/smalt sample/' task_smaltmap.pl
sed -i 's/$smalt map/smalt map/' task_smaltmap.pl
sed -i 's/$samtools sort/samtools sort/' task_smaltmap.pl
sed -i 's/$samtools view/samtools view/' task_smaltmap.pl
sed -i 's/$samtools rmdup/samtools rmdup/' task_smaltmap.pl
sed -i 's/$samtools reheader/samtools reheader/' task_smaltmap.pl
sed -i 's/$samtools index/samtools index/' task_smaltmap.pl

#install chrome
firefox https://www.google.com/chrome/browser/desktop/index.html
sudo dpkg -i google-chrome-stable_current_amd64.deb

#install SMRT link v4.0
# http://www.pacb.com/support/software-downloads/
export SMRT_USER="smrtanalysis" #| tee -a "${HOME}"/.bashrc
# source "${HOME}"/.bashrc
sudo addgroup smrtanalysis
sudo adduser --ingroup smrtanalysis smrtanalysis
sudo usermod -aG sudo smrtanalysis
#login as root
su -l $SMRT_USER
SMRT_USER="smrtanalysis"
export SMRT_ROOT=$PWD/pacbio/smrtlink
wget http://programs.pacificbiosciences.com/e/1652/lers-smrtlink-4-0-0-190159-zip/3rvmzd/507864561
unzip 507864561 #zip password: SmrT3chN
# mkdir -p /home/smrtanalysis/pacbio/data
# mkdir -p /home/smrtanalysis/pacbio/results
./smrtlink_4.0.0.190159.run --rootdir $SMRT_ROOT
#Start the services
$SMRT_ROOT/admin/bin/services-start
$SMRT_ROOT/admin/bin/services-status
#Validate the installation
$SMRT_ROOT/admin/bin/import-canneddata
$SMRT_ROOT/admin/bin/run-sat-services
google-chrome --enable-plugins http://localhost:9090/ #user/password: admin/admin  # do not change admin password.
#stop deamon
$SMRT_ROOT/admin/bin/services-stop



#install Pacbio SMRT analysis software
echo -e "#Pacbio SMRT analysis software"
export SMRT_ROOT="/opt/smrtanalysis"
export SMRT_USER="smrtanalysis"
export SMRT_GROUP="smrtanalysis" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
sudo addgroup smrtanalysis
sudo adduser --ingroup smrtanalysis smrtanalysis
sudo mkdir $SMRT_ROOT
sudo chown $SMRT_USER:$SMRT_GROUP $SMRT_ROOT
#login as root
su -l $SMRT_USER
wget http://files.pacb.com/software/smrtanalysis/2.3.0/smrtanalysis_2.3.0.140936.run
wget https://s3.amazonaws.com/files.pacb.com/software/smrtanalysis/2.3.0/smrtanalysis-patch_2.3.0.140936.p5.run
SMRT_ROOT=/opt/smrtanalysis
SMRT_USER=smrtanalysis
SMRT_GROUP=smrtanalysis
bash smrtanalysis_2.3.0.140936.run -p smrtanalysis-patch_2.3.0.140936.p5.run --rootdir $SMRT_ROOT
rm smrtanalysis*
#Start the SMRT Analysis daemons
$SMRT_ROOT/admin/bin/smrtportald-initd start
# $SMRT_ROOT/admin/bin/kodosd start
#Launch the SMRT Analysis portal
firefox http://localhost:8080/smrtportal &
#proceed to setup according to http://www.pacb.com/wp-content/uploads/SMRT-Analysis-Software-Installation-v2-3-0.pdf
#Stop the SMRT Analysis portal
$SMRT_ROOT/admin/bin/smrtportald-initd stop
#logout smrtanalysis user
exit

#Install SMRT view
cd "${HOME}"/prog
wget http://files.pacb.com/software/smrtanalysis/2.3.0/smrtview_2.3.0.140936_linux.tar.gz
tar zxvf smrtview_2.3.0.140936_linux.tar.gz
rm smrtview_2.3.0.140936_linux.tar.gz
cd smrtview_2.3.0/bin
./linux_configure

#install OpenMPI
sudo aptitude install libopenmpi-dev openmpi-bin openmpi-doc

#Install Ray
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/denovoassembler/Ray-2.3.1.tar.bz2
tar xvjf Ray-2.3.1.tar.bz2
rm Ray-2.3.1.tar.bz2
cd Ray-2.3.1
make PREFIX=/usr/local/bin
sudo make install
mpiexec -n $(nproc) Ray -o test -p test_1.fastq test_2.fastq -k 31

#Install SPAdes
cd "${HOME}"/prog
wget http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1.tar.gz
tar zxvf SPAdes-3.10.1.tar.gz
rm SPAdes-3.10.1.tar.gz
cd SPAdes-3.10.1
sudo PREFIX=/usr/local ./spades_compile.sh

#install sparsehash for ABySS
cd "${HOME}"/prog
git clone --recursive https://github.com/sparsehash/sparsehash.git
cd sparsehash
./configure
make
sudo make install

#Install ABySS
cd "${HOME}"/prog
git clone --recursive https://github.com/bcgsc/abyss.git
cd abyss
./autogen.sh
./configure --with-mpi=/usr/lib/openmpi --enable-maxk=128
make
sudo make install

# wget https://github.com/bcgsc/abyss/releases/download/1.9.0/abyss-1.9.0.tar.gz
# tar zxvf abyss-1.9.0.tar.gz
# rm abyss-1.9.0.tar.gz
# cd abyss-1.9.0
# ./configure
# make
# sudo make install

#install PEAR
cd "${HOME}"/prog
git clone --recursive https://github.com/xflouris/PEAR.git
cd PEAR
./autogen.sh
./configure
make
sudo make install

#install PhiSpy
sudo R -e 'install.packages("randomForest")'
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/phispy/phiSpyNov11_v2.3.zip
unzip phiSpyNov11_v2.3.zip
rm phiSpyNov11_v2.3.zip
rm -rf __MACOSX
cd phiSpyNov11_v2.3
make
chmod +x phiSpy.py
echo -e "\n#PhiSpy\nexport PATH=\$HOME/prog/phiSpyNov11_v2.3:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

# #install xgraph
# cd "${HOME}"/prog
# wget http://www.xgraph.org/linux/xgraph_4.30_linux64.tar.gz
# tar zxvf xgraph_4.30_linux64.tar.gz
# rm xgraph_4.30_linux64.tar.gz
# cd XGraph4.30_linux64
# sudo cp bin/xgraph /usr/local/bin
# echo -e "\n#XGraph\nexport PATH=\$HOME/prog/XGraph4.30_linux64:\$PATH" | tee -a "${HOME}"/.bashrc
# source "${HOME}"/.bashrc

# #install phage_finder
# sudo cpan 'Math::Round'
# sudo cpan 'PHAGE::Phage_subs_v2' #not found
# sudo cpan 'Getopt::Std'
# cd "${HOME}"/prog
# wget http://downloads.sourceforge.net/project/phage-finder/phage_finder_v2.1/phage_finder_v2.1.tar.gz
# tar zxvf phage_finder_v2.1.tar.gz
# rm phage_finder_v2.1.tar.gz
# cd phage_finder_v2.1


#install sublime text
# cd "${HOME}"/prog
# wget https://download.sublimetext.com/sublime-text_build-3114_amd64.deb
# sudo dpkg -i sublime-text_build-3114_amd64.deb
# rm sublime-text_build-3114_amd64.deb
wget -qO - https://download.sublimetext.com/sublimehq-pub.gpg | sudo apt-key add -
echo "deb https://download.sublimetext.com/ apt/stable/" | sudo tee /etc/apt/sources.list.d/sublime-text.list
sudo apt-get update
sudo apt-get install sublime-text

#install eclipse
sudo cpan 'File::Fetch'
cd "${HOME}"/prog
wget http://eclipse.mirror.rafal.ca/oomph/epp/mars/R2/eclipse-inst-linux64.tar.gz
tar zxvf eclipse-inst-linux64.tar.gz
rm eclipse-inst-linux64.tar.gz
cd eclipse-installer
./eclipse-inst
#change install fordel for /home/bioinfo/prog/eclipse

#Install EPIC in Eclipse
sudo cpan 'PadWalker'
http://www.epic-ide.org/updates/testing

#install PyDev in Eclipse
http://pydev.org/updates

#install shutter
sudo add-apt-repository ppa:shutter/ppa
sudo apt-get update
sudo apt-get install shutter

#install filezilla
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/filezilla/FileZilla_Client/3.18.0/FileZilla_3.18.0_x86_64-linux-gnu.tar.bz2
tar xvjf FileZilla_3.18.0_x86_64-linux-gnu.tar.bz2
rm FileZilla_3.18.0_x86_64-linux-gnu.tar.bz2
cd FileZilla3
echo -e "\n#PhiSpy\nexport PATH=\$HOME/prog/phiSpyNov11_v2.3:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install Display link adapter driver
sudo apt-get install dkms
cd "${HOME}"/prog
mkdir displaylink
cd displaylink
#download the driver from http://www.displaylink.com/downloads/ubuntu
unzip DisplayLink\ USB\ Graphics\ Software\ for\ Ubuntu\ 1.3.54.zip
rm DisplayLink\ USB\ Graphics\ Software\ for\ Ubuntu\ 1.3.54.zip
chmod +x displaylink-driver-1.3.54.run
sudo ./displaylink-driver-1.3.54.run
cd ..
rm -rf displaylink

#install LAST
cd "${HOME}"/prog
wget http://last.cbrc.jp/last-849.zip
unzip last-849.zip
rm last-849.zip
cd last-849
make -j
sudo make install

#install DIAMOND blast
cd "${HOME}"/prog
git clone --recursive https://github.com/bbuchfink/diamond.git
cd diamond
./build_simple.sh
sudo cp diamond /usr/local/bin/

#install MCMC decont
cd "${HOME}"/prog
git clone --recursive https://github.com/Lafond-LapalmeJ/MCSC_Decontamination.git
#uniref100 taxlist
cd /media/6tb_raid10/db
git clone https://github.com/GDKO/uniref_taxlist.git
cd uniref_taxlist
cat uniref100.taxlist.gz.part-a* | gunzip > uniref100.taxlist
ls -A | grep -F "part" | xargs rm -rf
cd ..
mv uniref_taxlist uniref
#uniref90
cd /media/6tb_raid10/db/uniref
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
pigz -d uniref90.fasta.gz
diamond makedb -v --threads $(nproc) --in uniref90.fasta -d uniref90
#ncbi taxonomy
cd "${HOME}"/prog/MCSC_Decontamination/data
ln -s /media/6tb_raid10/db/taxonomy/names.dmp names.dmp
ln -s /media/6tb_raid10/db/taxonomy/nodes.dmp nodes.dmp
# wget $(cat ncbi_taxo_ftp_link)
# tar zxvf taxdump.tar.gz
# rm taxdump.tar.gz

#install prinseq
sudo cpanm Getopt::Long Pod::Usage File::Temp Fcntl Cwd JSON Cairo \
    Statistics::PCA MIME::Base64 CGI File::Path IO::Uncompress::AnyUncompress \
    LWP::Simple File::Copy File::Basename
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz
tar zxvf prinseq-lite-0.20.4.tar.gz
rm prinseq-lite-0.20.4.tar.gz
cd prinseq-lite-0.20.4
chmod +x *.pl
sudo cp prinseq-graphs-noPCA.pl /usr/local/bin/prinseq-graphs-noPCA
sudo cp prinseq-graphs.pl /usr/local/bin/prinseq-graphs
sudo cp prinseq-lite.pl /usr/local/bin/prinseq-lite

#install seqprep
cd "${HOME}"/prog
git clone --recursive https://github.com/jstjohn/SeqPrep.git
cd SeqPrep
make
sudo cp SeqPrep /usr/local/bin

#install parallel
cd "${HOME}"/prog
wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
tar xvjf parallel-latest.tar.bz2
rm parallel-latest.tar.bz2
cd parallel*
./configure
make
sudo make install

#install blob tools
sudo apt-get install build-essential gfortran libatlas-base-dev python-pip python-dev python-matplotlib
sudo pip install --upgrade pip
sudo pip install virtualenv
cd "$HOME"
mkdir my_virtualenv
cd my_virtualenv
virtualenv blob_env
cd blob_env
source bin/activate
pip install matplotlib
pip install docopt
pip install ujson
cd "${HOME}"/prog
git clone --recursive https://github.com/DRL/blobtools.git
echo -e "\n#blobtools\nexport PATH=\$HOME/prog/blobtools:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
#test
cd blobtools
./blobtools create \
    -i test_files/assembly.fna \
    -c test_files/mapping_1.bam.cov \
    -t test_files/blast.out \
    -o test_files/test_output \
    --names /media/6tb_raid10/db/taxonomy/names.dmp \
    --nodes /media/6tb_raid10/db/taxonomy/nodes.dmp
./blobtools view \
    -i test_files/test_output.blobDB.json \
    -o test_files/
grep '^##' test_files/test_output.blobDB.table.txt
grep -v '^##' test_files/test_output.blobDB.table.txt | \
    column -t -s $'\t'
./blobtools blobplot \
  -i test_files/blobDB.json \
  -o test_files/ 
deactivate #to exit the virtual environment
# echo -e "\n#Blobtools\nexport PATH=\$HOME/prog/blobtools:\$PATH" | tee -a "${HOME}"/.bashrc

#install usearch 32bit (free)
#see http://www.drive5.com/usearch/download.html
cd "${HOME}"/prog
mkdir usearch
cd usearch
mv "${HOME}"/Downloads/usearch8.1.1861_i86linux32 ./usearch
chmod +x usearch
sudo cp usearch /usr/local/bin

#install PanPhlAn
#need usearch7 installed first.
#need bowtie2
cd "${HOME}"/prog
hg clone https://bitbucket.org/CibioCM/panphlan
echo -e "#\nPanPhlAn\nexport BOWTIE2_INDEXES=/media/6tb_raid10/db/panphlan\nexport PATH=\$HOME/prog/panphlan:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
#install python modules
sudo pip install sklearn
#install db
mkdir /media/6tb_raid10/db/panphlan
cd /media/6tb_raid10/db/panphlan
wget https://www.dropbox.com/sh/rfwb4i9m8s40iba/AAD8hB_JI1y5VIa8QLWMAReza/panphlan_senterica16.zip?dl=0
mv panphlan_senterica16.zip?dl=0 panphlan_senterica16.zip
unzip panphlan_senterica16.zip
rm panphlan_senterica16.zip
wget https://www.dropbox.com/sh/rfwb4i9m8s40iba/AAC4mIlAvNKKtC1x-5tIe93Aa/panphlan_ecoli16.zip?dl=0
mv panphlan_ecoli16.zip?dl=0 panphlan_ecoli16.zip
unzip panphlan_ecoli16.zip
rm panphlan_ecoli16.zip

#install MetaPhlAn2
sudo apt-get install mercurial
cd "${HOME}"/prog
hg clone https://bitbucket.org/biobakery/metaphlan2
echo -e "\n#MetaPhlAn2\nexport PATH=\$HOME/prog/metaphlan2:\$PATH\nexport mpa_dir=\$HOME/prog/metaphlan2" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc



#install MEGAN
cd "${HOME}"/prog
mkdir megan6
cd megan6
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_4_9.sh
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/prot_acc2tax-June2016.abin.zip
wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/nucl_acc2tax-June2016.abin.zip



#install FastTree
cd "${HOME}"/prog
mkdir FastTree
#see  http://www.microbesonline.org/fasttree/#Install
chmod +x FastTree*
sudo cp FastTreeMP /usr/local/bin
sudo cp FastTree /usr/local/bin


#install biom-format
cd "${HOME}"/prog
sudo pip install biom-format
sudo pip install h5py
sudo R
install.packages("biom")
quit()



#install STAMP
cd "${HOME}"/prog
sudo apt-get install libblas-dev liblapack-dev gfortran
sudo apt-get install freetype* python-pip python-dev python-numpy python-scipy python-matplotlib python-qt4
sudo pip install STAMP


#install Cytoscape
cd "${HOME}"/prog




#install ITSx
cd "${HOME}"/prog




#install PICRUSt
cd "${HOME}"/prog



#install RDPtools
cd "${HOME}"/prog
git clone --recursive https://github.com/rdpstaff/RDPTools.git
cd RDPTools
git clone --recursive https://github.com/rdpstaff/TaxonomyTree.git
cd ..
git submodule init
git submodule update
sudo make
echo -e "\n#RDPTools\nexport PATH=\$HOME/prog/RDPTools:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

# install RDP classifier
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.2.zip
unzip rdp_classifier_2.2.zip
rm rdp_classifier_2.2.zip
cd rdp_classifier_2.2
echo -e "\n#RDP classifier\nexport RDP_JAR_PATH=\$HOME/prog/rdp_classifier_2.2/rdp_classifier-2.2.jar" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install RDP paired-end
# cd "${HOME}"/prog
# wget http://rdp.cme.msu.edu/download/RDP_Assembler.tgz
# tar zxvf RDP_Assembler.tgz
# rm RDP_Assembler.tgz
# cd RDP_Assembler

#install Trimmomatic
cd "${HOME}"/prog
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
rm Trimmomatic-0.36.zip
cd Trimmomatic-0.36
#echo -e "\n#Trimmomatic\nexport PATH=\$HOME/prog/Trimmomatic-0.36/trimmomatic-0.36.jar:\$HOME/prog/Trimmomatic-0.36/adapters:\$PATH" | tee -a "${HOME}"/.bashrc
#source "${HOME}"/.bashrc

#install Flash
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz
tar zxvf FLASH-1.2.11.tar.gz
rm FLASH-1.2.11.tar.gz
cd FLASH-1.2.11
make
sudo cp flash /usr/local/bin


#install BLAT
cd "${HOME}"/prog
wget https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
unzip blatSrc35.zip
rm blatSrc35.zip
cd blatSrc
export MACHTYPE=x86_64
mkdir -p bin/"$MACHTYPE"
mkdir -p lib/"$MACHTYPE"
sed -i 's%${HOME}/bin/${MACHTYPE}%${HOME}/prog/blatSrc/bin/${MACHTYPE}%' inc/common.mk
make
sudo cp bin/"${MACHTYPE}"/* /usr/local/bin

#install LEfSe
sudo pip install rpy2
sudo R
install.packages("mvtnorm")
install.packages("modeltools")
install.packages("coin")
quit()_
cd "${HOME}"/prog
hg clone https://bitbucket.org/nsegata/lefse
#to be continued...

#install ViennaRNA
sudo apt-get install libperl-dev
cd "${HOME}"/prog
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_2_x/ViennaRNA-2.2.6.tar.gz
tar zxvf ViennaRNA-2.2.6.tar.gz
rm ViennaRNA-2.2.6.tar.gz
cd  ViennaRNA-2.2.6
./configure
make
sudo make install

#install Primer Prospector
sudo pip install --upgrade pip
sudo pip install cogent[all]
sudo pip install numpy
sudo pip install scipy
sudo pip install matplotlib
sudo pip install sphinx
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/pprospector/pprospector-1.0.1.tar.gz
tar zxvf pprospector-1.0.1.tar.gz
rm pprospector-1.0.1.tar.gz
cd pprospector-1.0.1
sudo python setup.py install


#install DiScRIBinATE
#see http://metagenomics.atc.tcs.com/binning/DiScRIBinATE/#download
cd "${HOME}"/prog
tar zxvf DiScRIBinATE.tar.gz
rm DiScRIBinATE.tar.gz
cd DiScRIBinATE
#to be continued...

#Download NCBI nr database
cd /media/6tb_raid10/db
mkdir nr
cd nr
for i in {00..52}; do
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr."${i}".tar.gz"
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr."${i}".tar.gz.md5"
done
md5sum -c *.md5 >> nr_md5check.txt
cat nr_md5check.txt
parallel --jobs $(nproc) tar zxf nr.{}.tar.gz ::: {00..52}
rm *.tar.gz*
#see ftp://ftp.ncbi.nlm.nih.gov/blast/db/
# parallel --jobs $(nproc) wget -nc "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.{}.tar.gz" ::: {00..52}
# parallel --jobs $(nproc) wget -nc "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.{}.tar.gz.md5" ::: {00..52}
# ls *.md5 | parallel --jobs $(nproc) md5sum -c *.md5 '>' nr_md5check.txt

#Download NCBI nt database
cd /media/6tb_raid10/db
mkdir nt
cd nt
for i in {00..38}; do
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt."${i}".tar.gz"
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt."${i}".tar.gz.md5"
done
md5sum -c *.md5 >> nt_md5check.txt
cat nt_md5check.txt
parallel --jobs $(nproc) tar zxf nt.{}.tar.gz ::: {00..38}
rm *.tar.gz*
# parallel --jobs $(nproc) wget -nc "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.{}.tar.gz" ::: {00..38}
# parallel --jobs $(nproc) wget -nc "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.{}.tar.gz.md5" ::: {00..38}

#install biopython
sudo pip install biopython
cd "${HOME}"/prog
git clone --recursive https://github.com/MG-RAST/DRISEE.git
echo -e "\n#DRISEE\nexport PATH=\$HOME/prog/DRISEE:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install muscle
cd "${HOME}"/prog
mkdir muscle_aligner
cd muscle
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar zxvf muscle3.8.31_i86linux64.tar.gz
rm muscle3.8.31_i86linux64.tar.gz
cd m

#install galaxy
cd "${HOME}"/prog
hg clone https://bitbucket.org/galaxy/galaxy-dist
cd galaxy-dist
./run.sh
#set admin account
#set toolshed account

#Upgrading all python packages with pip
sudo pip freeze --local | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U


#install docker
sudo apt-get update
sudo apt-get install \
    linux-image-extra-$(uname -r) \
    linux-image-extra-virtual
sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo apt-key fingerprint 0EBFCD88
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
sudo apt-get update
sudo apt-get install docker-ce
sudo docker run hello-world

#upgrade flash
sudo apt-get install flashplugin-installer

#install vcflib
cd "${HOME}"/prog
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
make -j
echo -e "\n#vcflib\nexport PATH=\$HOME/prog/vcflib/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install SnpSift
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip
# cd snpEff
# echo -e "\n#snpEff\nexport PATH=\$HOME/prog/snpEff/snpEff.jar:\$PATH" | tee -a "${HOME}"/.bashrc
# source "${HOME}"/.bashrc

#install Qualimap
sudo apt-get install libxml2-dev
sudo apt-get install libcurl4-openssl-dev
cd "${HOME}"/prog
wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.zip
unzip qualimap_v2.2.zip
rm qualimap_v2.2.zip
cd qualimap_v2.2
sudo Rscript scripts/installDependencies.r #### Error on dependency installation...
echo -e "\n#qualimap\nexport PATH=\$HOME/prog/qualimap_v2.2:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install Cap'n Proto
cd "${HOME}"/prog
git clone --recursive https://github.com/sandstorm-io/capnproto.git
cd capnproto/c++
autoreconf -i
./configure
make -j check
sudo make install

#install mash (MinHash)
# sudo apt-get install capnproto
sudo apt-get install libgsl-dev
cd "${HOME}"/prog
git clone --recursive https://github.com/marbl/Mash.git
cd Mash
./bootstrap.sh
./configure
make -j
sudo make install

#install Celera assembler
cd "${HOME}"/prog
wget http://downloads.sourceforge.net/project/wgs-assembler/wgs-assembler/wgs-8.3/wgs-8.3rc2-Linux_amd64.tar.bz2
tar xvjf wgs-8.3rc2-Linux_amd64.tar.bz2
rm wgs-8.3rc2-Linux_amd64.tar.bz2
echo -e "\n#Celera\nexport PATH=\$HOME/prog/wgs-8.3rc2/Linux-amd64/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install EDirect (ncbi)
cd "${HOME}"/prog
perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $ftp->login;
    $ftp->binary; $ftp->get("/entrez/entrezdirect/edirect.zip");'
unzip -u -q edirect.zip
rm edirect.zip
echo -e "\n#EDirect\nexport PATH=\$HOME/prog/edirect:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
cd edirect
sudo ./setup.sh

#install Ragout
cd "${HOME}"/prog
wget https://github.com/fenderglass/Ragout/releases/download/v1.2/ragout-1.2-linux-x86_64.tar.gz
tar zxvf ragout-1.2-linux-x86_64.tar.gz
rm ragout-1.2-linux-x86_64.tar.gz

#install CLARK
cd "${HOME}"/prog
wget http://clark.cs.ucr.edu/Download/CLARKV1.2.3.tar.gz
tar zxvf CLARKV1.2.3.tar.gz
rm CLARKV1.2.3.tar.gz
cd CLARKSCV1.2.3
./install.sh

#install mailutils
sudo apt-get install mailutils

#IGV install
sudo apt-get install igv

#install igvtools
cd "${HOME}"/prog
wget http://data.broadinstitute.org/igv/projects/downloads/igvtools_2.3.91.zip
unzip igvtools_2.3.91.zip
rm igvtools_2.3.91.zip

#install sra toolkit
cd "${HOME}"/prog
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz
tar zxvf sratoolkit.2.8.2-1-ubuntu64.tar.gz
rm sratoolkit.2.8.2-1-ubuntu64.tar.gz
echo -e "\n#sra toolkit\nexport PATH=\$HOME/prog/sratoolkit.2.8.2-1-ubuntu64/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc


#install parsnp
cd "${HOME}"/prog
wget https://github.com/marbl/parsnp/releases/download/v1.2/Parsnp-Linux64-v1.2.tar.gz
tar zxvf Parsnp-Linux64-v1.2.tar.gz
rm Parsnp-Linux64-v1.2.tar.gz
sudo cp Parsnp-Linux64-v1.2/parsnp /usr/local/bin

#install centrifuge
cd "${HOME}"/prog
git clone https://github.com/infphilo/centrifuge
cd centrifuge
make -j
sudo make install prefix=/usr/local
#install db
[ -d /media/6tb_raid10/db/centrifuge/taxonomy ] || mkdir -p /media/6tb_raid10/db/centrifuge/taxonomy
centrifuge-download -o /media/6tb_raid10/db/centrifuge/taxonomy taxonomy
[ -d /media/6tb_raid10/db/centrifuge/library ] || mkdir -p /media/6tb_raid10/db/centrifuge/library
centrifuge-download \
   -o /media/6tb_raid10/db/centrifuge/library \
   -d 'viral' \
   -m \
   -P $(nproc) \
   refseq > /media/6tb_raid10/db/centrifuge/seqid2taxid.map
cat /media/6tb_raid10/db/centrifuge/library/*/*.fna > /media/6tb_raid10/db/centrifuge/input-sequences.fna
today=$(date +%Y-%m-%d)
centrifuge-build \
    -p $(nproc) \
    --conversion-table /media/6tb_raid10/db/centrifuge/seqid2taxid.map \
    --taxonomy-tree /media/6tb_raid10/db/centrifuge/taxonomy/nodes.dmp \
    --name-table /media/6tb_raid10/db/centrifuge/taxonomy/names.dmp \
    /media/6tb_raid10/db/centrifuge/input-sequences.fna \
    /media/6tb_raid10/db/centrifuge/"${today}"_refseq_viral_complete_genome
#Download nt for centrifuge ready to go
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/nt.tar.gz


##################### Roary ######################

#Install MCL
cd "${HOME}"/prog
wget http://www.micans.org/mcl/src/mcl-14-137.tar.gz
tar zxvf mcl-14-137.tar.gz 
rm mcl-14-137.tar.gz
cd mcl-14-137
./configure --enable-blast
make -j
sudo make install

#install Roary
sudo apt-get install cpanminus
sudo cpanm Array::Utils Bio::Perl Exception::Class File::Basename File::Copy File::Find::Rule File::Grep File::Path File::Slurper File::Spec File::Temp File::Which FindBin Getopt::Long Graph Graph::Writer::Dot List::Util Log::Log4perl Moose Moose::Role Text::CSV PerlIO::utf8_strict 
# sudo cpanm -f Bio::Roary
#roary
cd "${HOME}"/prog
git clone --recursive https://github.com/sanger-pathogens/Roary.git
echo -e "\n#Roary\nexport PATH=\$HOME/prog/Roary/bin:\$PATH\nexport PERL5LIB=\$PERL5LIB:\$HOME/prog/Roary/lib" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install scoary
sudo pip install ete3 six
sudo pip install scoary

###################################################

# Install Gubbins
# sudo apt-get install gubbins
cd "${HOME}"/prog
git clone --recursive https://github.com/sanger-pathogens/gubbins.git
cd gubbins
sudo apt-get install python3-dev
autoreconf -i
./configure
make -j
sudo make install
cd python
sudo python3 setup.py install

#install sabre
cd "${HOME}"/prog
git clone --recursive https://github.com/najoshi/sabre.git
cd sabre
make -j
echo -e "\n#sabre\nexport PATH=\$HOME/prog/sabre:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install poretools
sudo apt-get install git python-setuptools python-dev cython libhdf5-serial-dev  # dependencies
cd "${HOME}"/prog
git clone https://github.com/arq5x/poretools
cd poretools
sudo python setup.py install

#install nanopolish
cd "${HOME}"/prog
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish
make -j
echo -e "\n#nanopolish\nexport PATH=\$HOME/prog/nanopolish:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install racon
cd "${HOME}"/prog
git clone https://github.com/isovic/racon.git
cd racon
make modules
make tools
make -j
sudo cp bin/racon /usr/local/bin/
sudo cp tools/minimap/minimap /usr/local/bin/

# install sbt
echo "deb https://dl.bintray.com/sbt/debian /" | sudo tee -a /etc/apt/sources.list.d/sbt.list
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823
sudo apt-get update
sudo apt-get install sbt

#install Pilon
cd "${HOME}"/prog
git clone --recursive https://github.com/broadinstitute/pilon.git
cd pilon
./build.sh
#fix the ln line

#Install Jaspa
cd "${HOME}"/prog
wget https://github.com/mdcao/japsa/releases/download/v1.7-02a/JapsaRelease.tar.gz
tar zxvf JapsaRelease.tar.gz
cd JapsaRelease
sudo bash install.sh

#Install npScaf
cd "${HOME}"/prog
git clone https://github.com/mdcao/japsa
cd japsa
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
sudo make install \
   INSTALL_DIR=/usr/local/bin \
   MXMEM="${mem}"g

#install phage_typing
cd "${HOME}"/prog
git clone https://github.com/OLF-Bioinformatics/phage_typing
echo -e "\n#phage_typing\nexport PATH=\$HOME/prog/phage_typing:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install snp_analysis (USDA)
cd "${HOME}"/prog
git clone https://github.com/OLF-Bioinformatics/snp_analysis
echo -e "\n#snp_analysis\nexport PATH=\$HOME/prog/snp_analysis:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install MAFFT
cd "${HOME}"/prog
wget http://mafft.cbrc.jp/alignment/software/mafft-7.310-with-extensions-src.tgz
tar zxvf mafft-7.310-with-extensions-src.tgz
rm mafft-7.310-with-extensions-src.tgz
mv mafft-7.310-with-extensions mafft
cd mafft/core
make clean
make -j
sudo make install

#install Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
#follow on-screen instruction
rm Miniconda3-latest-Linux-x86_64.sh

#install SISTR
sudo pip install --upgrade pip
conda create -n sistr python=3.4.5
source activate sistr
pip install sistr_cmd
source activate

#install SeqSero
cd "${HOME}"/prog
git clone --recursive https://github.com/denglab/SeqSero.git
echo -e "\n#SeqSero\nexport PATH=\$HOME/prog/SeqSero:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc


#install gff utilities
cd "${HOME}"/prog
git clone --recursive https://github.com/gpertea/gclib.git
git clone --recursive https://github.com/gpertea/gffread.git
cd gffread
make
sudo cp gffread /usr/local/bin

#install ncbi-tools (sequin)
sudo apt install ncbi-tools-x11

#install OPERA
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/operasf/OPERA-LG%20version%202.0.6/OPERA-LG_v2.0.6.tar.gz
tar zxvf OPERA-LG_v2.0.6.tar.gz
rm OPERA-LG_v2.0.6.tar.gz
cd OPERA-LG_v2.0.6
make install
echo -e "\n#OPERA\nexport PATH=\$HOME/prog/OPERA-LG_v2.0.6/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install BESST
sudo apt-get install python-networkx python-pysam
cd "${HOME}"/prog
git clone https://github.com/ksahlin/BESST.git
cd BESST
sudo python setup.py install
# sudo pip install BESST
# runBESST

#install DISCOVAR de novo
cd "${HOME}"/prog
wget https://github.com/jemalloc/jemalloc/releases/download/4.5.0/jemalloc-4.5.0.tar.bz2
tar jxvf jemalloc-4.5.0.tar.bz2
rm jemalloc-4.5.0.tar.bz2
cd jemalloc-4.5.0
./configure
make -j
sudo make install
cd "${HOME}"/prog
wget ftp://ftp.broadinstitute.org/pub/crd/DiscovarDeNovo/latest_source_code/LATEST_VERSION.tar.gz
tar zxvf LATEST_VERSION.tar.gz
rm LATEST_VERSION.tar.gz
cd discovardenovo-52488
./configure
make -j
sudo make install

#install SSPACE
cd "${HOME}"/prog
wget https://www.baseclear.com/base/download/41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
tar zxvf 41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz

#install pyScaf
sudo pip install -U FastaIndex
cd "${HOME}"/prog
wget http://last.cbrc.jp/last-852.zip
unzip last-852.zip
rm last-852.zip
cd last-852
make -j
sudo make install
cd "${HOME}"/prog
git clone https://github.com/lpryszcz/pyScaf.git

#install pycharm
sudo cpanm Devel::Camelcadedb
cd "${HOME}"/prog
wget https://download.jetbrains.com/python/pycharm-community-2017.1.2.tar.gz
tar zxvf pycharm-community-2017.1.2.tar.gz
rm pycharm-community-2017.1.2.tar.gz
cd pycharm-community-2017.1.2/bin
./pycharm.sh


#install mummer
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
tar zxvf MUMmer3.23.tar.gz
rm MUMmer3.23.tar.gz
cd MUMmer3.23
make check 
sudo make install
echo -e "\n#MUMmer3\nexport PATH=\$HOME/prog/MUMmer3.23:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc


#install Metassembler - requires mate pairs
sudo cpan Statistics::Descriptive
sudo pip install argparse
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/metassembler/v1.5/Metassembler.1.5.tar.gz
tar zxvf Metassembler.1.5.tar.gz
rm Metassembler.1.5.tar.gz
cd Metassembler
make install
echo -e "\n#Metassembler\nexport PATH=\$HOME/prog/Metassembler/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install CAMSA
sudo pip install camsa


#install blat
cd "${HOME}"/prog
wget http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
unzip blatSrc35.zip
rm blatSrc35.zip
cd blatSrc
MACHTYPE=$(echo $MACHTYPE | cut -d "-" -f 1)
export MACHTYPE
mkdir -p lib/$MACHTYPE
sed -i 's%BINDIR = ${HOME}/bin/${MACHTYPE}%BINDIR = /usr/local/bin%' ./inc/common.mk
sudo make


#install BLAT
cd "${HOME}"/prog
wget http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
unzip blatSrc35.zip
rm blatSrc35.zip
cd blatSrc
MACHTYPE=$(echo $MACHTYPE | cut -d "-" -f 1)
export MACHTYPE
mkdir -p lib/$MACHTYPE
sed -i 's%BINDIR = ${HOME}/bin/${MACHTYPE}%BINDIR = /usr/local/bin%' ./inc/common.mk
sudo make
sudo cpan DBI
#install Qt4
cd "${HOME}"/prog
wget http://download.qt.io/official_releases/qt/5.8/5.8.0/single/qt-everywhere-opensource-src-5.8.0.tar.gz
tar zxvf qt-everywhere-opensource-src-5.8.0.tar.gz
rm qt-everywhere-opensource-src-5.8.0.tar.gz
cd qt-everywhere-opensource-src-5.8.0
./configure  # this is interactive
make -j
sudo make install

# wget http://download.qt.io/official_releases/online_installers/qt-unified-linux-x64-online.run
# chmod +x qt-unified-linux-x64-online.run
# sudo ./qt-unified-linux-x64-online.run


#install AMOS
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/amos/amos/3.1.0/amos-3.1.0.tar.gz
tar zxvf amos-3.1.0.tar.gz
rm amos-3.1.0.tar.gz
cd amos-3.1.0
./configure
make -j
sudo make install

#install circlator
echo -e "\n#Circlator\nexport CIRCLATOR_SPADES=/usr/local/bin/spades.py"  | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
#install pip3
conda install circlator

#install canu
cd "${HOME}"/prog
git clone https://github.com/marbl/canu.git
cd canu/src
make -j
echo -e "\n#Canu\nexport PATH=\$HOME/prog/canu/Linux-amd64/bin:\$PATH" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

# #install garm
# sudo cpan Parallel::ForkManager List::MoreUtils
# cd "${HOME}"/prog
# wget https://downloads.sourceforge.net/project/garm-meta-assem/GARM_v0.7.5.tar.gz
# tar zxvf GARM_v0.7.5.tar.gz
# rm GARM_v0.7.5.tar.gz
# cd GARM_v0.7.5
# echo -e "
# # GARM
# export GARMBIN=\$HOME/prog/GARM_v0.7.5/bin
# export GARMLIB=\$HOME/prog/GARM_v0.7.5/lib
# export MUMBIN=\$HOME/prog/GARM_v0.7.5/MUMmer3.22
# export AMOSBIN=\$HOME/prog/GARM_v0.7.5/amos-3.0.0/bin
# export AMOSLIB=\$HOME/prog/GARM_v0.7.5/amos-3.0.0/lib
# export PATH=\$HOME/prog/GARM_v0.7.5:\$PATH" | tee -a "${HOME}"/.bashrc
# source "${HOME}"/.bashrc

#install zorro
# http://lge.ibi.unicamp.br/zorro/
cd "${HOME}"/prog
wget http://lge.ibi.unicamp.br/zorro/downloads/Zorro2.2.tar.gz
tar zxvf Zorro2.2.tar.gz
rm Zorro2.2.tar.gz
cd Zorro2.2

#install CISA
# Requires MUMmer and blast+
cd "${HOME}"/prog
wget http://sb.nhri.org.tw/CISA/upload/en/2014/3/CISA_20140304-05194132.tar
tar xvf CISA_20140304-05194132.tar
rm CISA_20140304-05194132.tar
cd CISA

#install augustus
# http://bioinf.uni-greifswald.de/augustus/
sudo apt-get install libboost-iostreams-dev zlib1g-dev libgsl-dev libboost-graph-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev
cd "${HOME}"/prog
wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.3.tar.gz
tar zxvf augustus-3.2.3.tar.gz
rm augustus-3.2.3.tar.gz
cd augustus-3.2.3
echo -e "\n#Augustus\nexport PATH=\$HOME/prog/augustus-3.2.3/bin:\$HOME/prog/augustus-3.2.3/scripts:\$PATH"  | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install GeneMark-ES
sudo cpan YAML Hash::Merge Logger::Simple Parallel::ForkManager
cd "${HOME}"/prog
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_kPP34/gm_et_linux_64.tar.gz
tar zxvf gm_et_linux_64.tar.gz
rm gm_et_linux_64.tar.gz
cd 

#install ITSx
cd "${HOME}"/prog
wget http://microbiology.se/sw/ITSx_1.0.11.tar.gz
tar zxvf ITSx_1.0.11.tar.gz
rm ITSx_1.0.11.tar.gz
cd ITSx_1.0.11/ITSx_db/HMMs
rm *.hmm.*
find . -maxdepth 1 -type f | parallel 'hmmpress {}'

#install EMBOSS
cd "${HOME}"/prog
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
tar zxvf EMBOSS-6.6.0.tar.gz
rm EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0
./configure
make -j
sudo ldconfig
sudo make install

#install bedtools
cd "${HOME}"/prog
git clone --recursive https://github.com/arq5x/bedtools2.git
cd bedtools2
make -j
sudo cp bin/* /usr/local/bin

#install blast2go
cd "${HOME}"/prog
wget https://www.blast2go.com/downloads/software/blast2go/latest/4_1/Blast2GO_unix_4_1_x64.zip
unzip Blast2GO_unix_4_1_x64.zip
sudo bash Blast2GO_unix_4_1_x64.sh
rm Blast2GO_unix_4_1_x64.*


#install interproscan
# https://github.com/ebi-pf-team/interproscan/wiki
cd "${HOME}"/prog
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.24-63.0/interproscan-5.24-63.0-64-bit.tar.gz
tar pzxvf interproscan-5.24-63.0-64-bit.tar.gz
rm interproscan-5.24-63.0-64-bit.tar.gz
cd interproscan-5.24-63.0-64-bit/data
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-11.1.tar.gz
tar pzxvf panther-data-11.1.tar.gz
rm panther-data-11.1.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/lookup_service/lookup_service_5.24-63.0.tar.gz
which signalp
ln -s /home/bioinfo/prog/signalp-4.1/signalp bin/signalp/4.1/signalp
#install phobius
cd "${HOME}"/prog
mkdir phobius
cd phobius
wget http://phobius.binf.ku.dk/download/homologhmm_1.05.tar.gz
tar zxvf homologhmm_1.05.tar.gz
rm homologhmm_1.05.tar.gz
cp * ../interproscan-5.24-63.0/bin/phobius/1.01/
mv phobius.pl
wget https://github.com/ElofssonLab/TOPCONS2/blob/master/topcons2_webserver/predictors/PolyPhobius/phobius.model
which decodeanhmm.Linux_x86_64
bin/tmhmm/2.0c/decodeanhmm, data/tmhmm/2.0c/TMHMM2.0c.model
ln -s /usr/local/bin/decodeanhmm.Linux_x86_64 decodeanhmm

ln -s /home/bioinfo/prog/tmhmm-2.0c/lib/TMHMM2.0.model data/tmhmm/2.0c/TMHMM2.0c.model



#install local blast2go
# https://www.blast2go.com/start-blast2go-2-8
# https://www.blast2go.com/support/blog/22-blast2goblog/110-local-blast2go-database-installation
# http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/install/blast2godb.htm
sudo apt-get install mysql-server
cd "${HOME}"/prog
mkdir blast2go_local
cd blast2go_local
wget https://www.blast2go.com/webstart/blast2go10000.jnlp
cd local_b2g_db
mkdir -p /media/6tb_raid10/db/blast2go
cd /media/6tb_raid10/db/blast2go
echo -e "\
https://www.blast2go.com/images/b2g_support/local_b2g_db.zip
http://archive.geneontology.org/latest-full/go_monthly-assocdb-data.gz
ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz
ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz\
" > download_links.txt
cat download_links.txt | parallel 'wget {}'
unzip *.zip
pigz -d *.gz
sudo nano /etc/mysql/mysql.conf.d/mysqld.cnf
# below [mysqld] add: sql-mode="NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION"
sudo service mysql restart
mysql -h localhost -u root -p < local_b2g_db/b2gdb.sql
mysql -h localhost -u root -p -e "GRANT ALL ON b2gdb.* TO 'blast2go'@'localhost' IDENTIFIED BY 'blast4it';"
mysql -h localhost -u root -p -e "FLUSH PRIVILEGES;"
mysql -s -h localhost -u root -p b2gdb < go_monthly-assocdb-data
mysql -h localhost -u root -p b2gdb -e "LOAD DATA LOCAL INFILE '/media/6tb_raid10/db/blast2go/gene2accession' INTO TABLE gene2accession FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';"
mysql -h localhost -u root -p b2gdb -e "LOAD DATA LOCAL INFILE '/media/6tb_raid10/db/blast2go/gene_info' INTO TABLE gene_info FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';"
cd local_b2g_db
java -cp .:mysql-connector-java-5.0.8-bin.jar: ImportIdMapping /media/6tb_raid10/db/blast2go/idmapping.tb localhost b2gdb blast2go blast4it

#install blast (not blast+)
cd "${HOME}"/prog
mkdir blast
cd blast
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.9/blast-2.2.9-amd64-linux.tar.gz
tar zxvf blast-2.2.9-amd64-linux.tar.gz
rm blast-2.2.9-amd64-linux.tar.gz
echo -e "\n#BLAST-old\nexport PATH=\$HOME/prog/blast:\$PATH"  | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install PSORTb
# cd "${HOME}"/prog
# wget http://www.psort.org/download/libpsortb-1.0.tar.gz
# tar zxvf libpsortb-1.0.tar.gz
# rm libpsortb-1.0.tar.gz
# cd libpsortb-1.0
# ./configure
# make -j
# sudo make install
# sudo ldconfig
# sudo cpan XML::RPC
# #install blast (not blast+)
# cd "${HOME}"/prog
# wget http://www.psort.org/download/bio-tools-psort-all.3.0.4.tar.gz
# tar zxvf  bio-tools-psort-all.3.0.4.tar.gz
# rm bio-tools-psort-all.3.0.4.tar.gz
# cd bio-tools-psort-all
# mkdir blast
# cd blast
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.9/blast-2.2.9-ia64-linux.tar.gz
# tar zxvf blast-2.2.9-ia64-linux.tar.gz
# rm blast-2.2.9-ia64-linux.tar.gz
# cd ..
# mkdir pftools
# cd pftools
# wget ftp://ftp.lausanne.isb-sib.ch/pub/software/unix/pftools/pft2.3/executables/linux_x86_elf/rh7x/*
# chmod +x *
# cd ..
# cd algorithm-hmm
# perl Makefile.PL

# make
# make test
# make install
# #make blast db

# #add to PATH
# echo -e "\n#export PSORT_ROOT=\$HOME/" | tee -a "${HOME}"/.bashrc
# echo -e "export PSORT_PFTOOLS=" | tee -a "${HOME}"/.bashrc
# echo -e "export BLASTDIR=" | tee -a "${HOME}"/.bashrc
# source "${HOME}"/.bashrc


#install Phobius
cd "${HOME}"/prog
mkdir phobius
cd phobius
wget http://phobius.sbc.su.se/download/homologhmm_1.05.tar.gz
tar zxvf homologhmm_1.05.tar.gz
rm homologhmm_1.05.tar.gz
chmod +x jphobius blastget homologhmm
# get the model file by registering:
# http://software.sbc.su.se/cgi-bin/request.cgi?project=phobius


#install RepeatMasker
# use with hmmer
cd "${HOME}"/prog
sudo cpan Text::Soundex
# see http://tandem.bu.edu/trf/trf409.linux64.download.html
mv trf409.linux64 trf
sudo mv trf /usr/local/bin
wget http://www.repeatmasker.org/RepeatMasker-open-4-0-7.tar.gz
tar zxvf RepeatMasker-open-4-0-7.tar.gz
rm RepeatMasker-open-4-0-7.tar.gz
cd RepeatMasker/Librairies
wget http://www.dfam.org/web_download/Current_Release/Dfam.hmm.gz
pigz -d -f Dfam.hmm.gz
cd ..
# optional - Add RepBase RepeatMasker Edition database
wget RepBaseRepeatMaskerEdition-########.tar.gz
tar zxvf RepBaseRepeatMaskerEdition-########.tar.gz
rm RepBaseRepeatMaskerEdition-########.tar.gz
perl ./configure
echo -e "\n#RepeatMasker\nexport PATH=\$HOME/prog/RepeatMasker:\$PATH"  | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc

#install miropeats
#install ICAass first
cd "${HOME}"/prog
wget http://www.littlest.co.uk/software/pub/bioinf/freeold/icatools_2.5.tgz
tar zxvf icatools_2.5.tgz
rm icatools_2.5.tgz
cd icatools_2.5
make miropeats
echo -e "\n#ICAass\nexport PATH=\$HOME/prog/icatools_2.5/bin:\$PATH"  | tee -a "${HOME}"/.bashrc
#miropeats
source "${HOME}"/.bashrc
cd "${HOME}"/prog
wget http://www.littlest.co.uk/software/pub/bioinf/freeold/miropeats_2.02.tgz
tar zxvf miropeats_2.02.tgz
rm miropeats_2.02.tgz
cd miropeats_2.02
# ./miropeats -o ~/Desktop/Lw_mito_repeats.ps -USA -seq /home/bioinfo/analyses/Lw_assembly/mito/Lw_mito.fna > ~/Desktop/Lw_mito_repeats.txt


#install Pollux
cd "${HOME}"/prog
git clone --recursive https://github.com/emarinier/pollux.git
cd pollux
make -j
sudo mv pollux /usr/local/bin

#install Lighter
cd "${HOME}"/prog
git clone --recursive https://github.com/mourisl/Lighter.git
cd Lighter
make -j
sudo mv lighter /usr/local/bin


#install asn2gb
cd "${HOME}"/prog
wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/asn2gb/linux64.asn2gb.gz
pigz -d linux64.asn2gb.gz
mv linux64.asn2gb asn2gb
chmod +x asn2gb
sudo mv asn2gb /usr/local/bin

#install abricate
cd "${HOME}"/prog
git clone --recursive https://github.com/tseemann/abricate.git
cd abricate

#install circulator
cd "${HOME}"/prog
git clone --recursive https://github.com/sanger-pathogens/circlator.git
cd circlator

#install Phylip
# http://evolution.genetics.washington.edu/phylip/
# v3.696
cd "${HOME}"/prog
wget http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz
tar zxvf phylip-3.696.tar.gz
rm phylip-3.696.tar.gz
cd phylip-3.696/src
mv Makefile.unx Makefile
sudo make install

#install PlasmidFinder
cd "${HOME}"/prog
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
cd plasmidfinder
./INSTALL_DB database
./UPDATE_DB database
./VALIDATE_DB database

#install Aspera plugin for firefox
conda install libgcc
# http://downloads.asperasoft.com/connect2//
sudo cpan Try::Tiny::Retry

#install KAT
# conda config --add channels defaults
# conda config --add channels conda-forge
# conda config --add channels bioconda
# conda create --name kat python=3 kat
# source activate kat
# conda install -n kat gnuplot
# source deactivate

sudo apt-get install libboost-all-dev
cd "${HOME}"/prog
git clone --recursive https://github.com/TGAC/KAT.git
cd KAT
./autogen.sh
./configure
make -j
make check
sudo make install

#install Fast5-to-Fastq
cd "${HOME}"/prog
sudo pip3 install h5py
git clone https://github.com/rrwick/Fast5-to-Fastq
#Fast5-to-Fastq/fast5_to_fastq.py --help

#install Porechop
cd "${HOME}"/prog
git clone https://github.com/rrwick/Porechop.git
cd Porechop
sudo python3 setup.py install

#install tablet
cd "${HOME}"/prog
wget http://bioinf.hutton.ac.uk/tablet/installers/tablet_linux_x64_1_16_09_06.sh
chmod +x tablet_linux_x64_1_16_09_06.sh
sudo bash tablet_linux_x64_1_16_09_06.sh
rm tablet_linux_x64_1_16_09_06.sh


#installing poRe Parallel GUI
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/rpore/0.24/poRe_0.24.tar.gz

#install HPG-pore (only for old 2D nanopore data)
sudo apt-get install maven
cd "${HOME}"/prog
git clone https://github.com/opencb/hpg-pore.git
cd hpg-pore
./build.sh
echo -e "\n#HPG-pore\nexport PATH=\$HOME/prog/hpg-pore:\$PATH" | tee -a "${HOME}"/.bashrc
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/prog/hpg-pore" | tee -a "${HOME}"/.bashrc
source "${HOME}"/.bashrc
#./hpg-pore.sh -h

#install nanonet
apt-get install python-pyopencl
cd "${HOME}"/prog
git clone --recursive https://github.com/nanoporetech/nanonet.git
cd nanonet
sudo python setup.py install

#install genomescope
#Need jellyfish
cd "${HOME}"/prog
git clone --recursive https://github.com/schatzlab/genomescope.git

#install LmCGST
cd "${HOME}"/prog
mkdir LmCGST
cd LmCGST
wget https://downloads.sourceforge.net/project/lmcgst/LmCGSTv1.0.tar.gz
tar zxvf LmCGSTv1.0.tar.gz
rm LmCGSTv1.0.tar.gz
#ask Nick Petronella for v2.0 Perl script and new README. Database are still good

#install LoRMA (PacBio read self-correction)
cd "${HOME}"/prog
wget https://www.cs.helsinki.fi/u/lmsalmel/LoRMA/LoRMA-0.4.tar.gz
tar zxvf LoRMA-0.4.tar.gz
rm LoRMA-0.4.tar.gz
cd LoRMA-0.4
mkdir build;  cd build;  cmake ..;  make

#install jabba (hybrid reads correction)
cd "${HOME}"/prog
git clone https://github.com/aminallam/karect.git
cd karect
make -j
sudo cp karect /usr/local/bin

cd "${HOME}"/prog
git clone https://github.com/sparsehash/sparsehash.git
cd sparsehash
./configure
make -j
sudo make install

cd "${HOME}"/prog
git clone https://github.com/biointec/brownie.git
cd brownie
mkdir build
mkdir bin
cd build
cmake -DMAXKMERLENGTH=75 ..
make brownie
sudo cp src/brownie /usr/local/bin
# cp src/brownie ../bin/brownie

cd "${HOME}"/prog
git clone --recursive https://github.com/biointec/jabba.git
cd jabba
mkdir -p build
cd build
cmake ../
make
mkdir -p ../bin
sudo cp src/jabba /usr/local/bin
# cp -b src/jabba ../bin/jabba
# echo -e "\n#jabba\nexport PATH=\$HOME/prog/jabba/bin:\$PATH" | tee -a "${HOME}"/.bashrc
# source "${HOME}"/.bashrc


#install poligraph dependencies
#install cortex_var
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/cortexassembler/cortex_var/latest/CORTEX_release_v1.0.5.21.tgz
#install vcftools
cd "${HOME}"/prog
wget https://downloads.sourceforge.net/project/vcftools/vcftools_0.1.13.tar.gz
#install stampy
cd "${HOME}"/prog
#install poligraph
cd "${HOME}"/prog
git clone https://github.com/iqbal-lab/poligraph




#install Neptune
