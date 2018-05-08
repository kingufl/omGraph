#!/bin/bash
#set -xe
mkdir -p kmc_temp

D=100000
RES_ENZ="CGGACCG"

KMER_SIZE=31
READ_FILE="/ufrc/boucher/kingdgp/cosmo/cosmo/experiments/ecoli/corrected/pe_ecoli.fastq"
FA_FQ="fq"

echo $READ_FILE

echo /ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc -ci2 -$FA_FQ -k$KMER_SIZE -cs300 $READ_FILE $KMER_SIZE.kmc kmc_temp

/ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc -ci2 -$FA_FQ -k$KMER_SIZE -cs300 $READ_FILE $KMER_SIZE.kmc kmc_temp

echo /ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc_tools sort $KMER_SIZE.kmc $KMER_SIZE.kmc.sorted 

/ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc_tools sort $KMER_SIZE.kmc $KMER_SIZE.kmc.sorted 

echo /ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc_dump $KMER_SIZE.kmc.sorted $KMER_SIZE.kmc.sorted.dump 

/ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc_dump $KMER_SIZE.kmc.sorted $KMER_SIZE.kmc.sorted.dump 

echo "$KMER_SIZE.kmc.sorted" > $KMER_SIZE.reads.list

numactl --interleave=all /usr/bin/time -v /ufrc/boucher/kingdgp/cosmo/vari/cosmo/cosmo-build -d $KMER_SIZE.reads.list

cd ..

./vari_rest_kmers-o $RES_ENZ1 /ufrc/boucher/kingdgp/cosmo/cosmo/experiments/$KMER_SIZE.list.packed

cat rest_* > restriction_nodes$KMER_SIZE

exit 1

mkdir edges$KMER_SIZE

./dp_ext 100000 $KMER_SIZE

exit 1

ls -1 --color=no *.fastq |xargs -l -i echo "/ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc -ci2 -fq -k46 -cs300 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.fastq |xargs -l -i echo "/ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh

ls -1 --color=no *.fastq |xargs -l -i echo "/ufrc/boucher/kingdgp/cosmo/cosmo/3rd_party_src/KMC/bin/kmc_dump {}.kmc.sorted {}.kmc.sorted.dump " >kmerdumpsort.sh
source kmerdumpsort.sh

ls -1 --color=no *.fastq |xargs -l -i echo "{}.kmc.sorted" > ecoli6_kmc2_list

#ls -1 --color=no *.fastq |xargs -l -i echo "./find_restriction_nodes {}.kmc.sorted.dump" > find_res_nodes.sh
#source find_res_nodes.sh 

numactl --interleave=all /usr/bin/time -v /ufrc/boucher/kingdgp/cosmo/cosmo/cosmo-pack -k ecoli6_kmc2_list
