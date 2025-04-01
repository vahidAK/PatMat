#!/bin/bash -l
# 1- Installing patmat and plus some additional tools such as longphase, whatshap, samtools, sniffles, bowtie2
###################################################################################################################
git clone https://github.com/vahidAK/PatMat.git && cd PatMat &&\
chmod -x */* && chmod +x */{patmat.py,strandseq_phase.R} && \
mamba env create -f patmat_env.yml && \
conda activate patmat && \
env_path=$(mamba info | grep -i 'active env location' | cut -d'/' -f2- | awk '{print "/"$0"/bin"}') && \
r_libPaths=$(mamba info | grep -i 'active env location' | cut -d'/' -f2- | awk '{print "/"$0"/lib/R/library"}') && \
ln -s $PWD/{patmat,Strand-seq}/* $env_path && \
mkdir third_parties && cd third_parties && \
wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v3.0.0/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu.tar.gz && \
tar -xzf pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu.tar.gz && rm *.gz && \
chmod +x pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores && \
ln -s $PWD/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores $env_path
patmat.py --version
# 2- Installing ashleys-qc
###################################################################################################################
wget https://github.com/friendsofstrandseq/ashleys-qc/archive/refs/tags/v0.2.0.tar.gz && \
tar -xzf v0.2.0.tar.gz && rm v0.2.0.tar.gz && cd ashleys-qc-0.2.0 && \
sed -i 's/name: ashleys/name: ashleys_patmat-wf/g' environment/ashleys_env.yml && \
mamba env create -f environment/ashleys_env.yml && \
conda activate ashleys_patmat-wf && \
python setup.py install && \
chmod +x bin/ashleys.py && cd ../ \
ashleys.py -h

# 3- Installing clair3. Note that you need to download clair3 models yourself
################################################################################################################### 
# make sure channels are added in mamba
mamba config --add channels defaults
mamba config --add channels bioconda
mamba config --add channels mamba-forge

# create mamba environment named "clair3"
mamba create -n clair3_patmat-wf -c bioconda clair3=1.0.10 python=3.9.0 -y
conda activate clair3_patmat-wf
mkdir clair3_models && cd clair3_models && \
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz && \
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r1041_e82_400bps_sup_v420.tar.gz && \
tar -xzf clair3_models.tar.gz && tar -xzf r1041_e82_400bps_sup_v420.tar.gz && rm *.gz
run_clair3.sh --version

