## omGraph

omGraph is a software tool for aligning optical maps to de-bruijn graphs built from sequence reads. It requires the VARI software for building the de-buijn graph. The steps required for building the Vari software are mentioned below.

## Building notes

Five third party packages are required for VARI. All should be cloned within the 3rd_party_src directory.  Any 3rd party software may change in incompatible ways.  The revisions known to work for the published results are included.


1. KMC2 --  'git clone https://github.com/refresh-bio/KMC' (commit f090276855a3f7c0b14e9f3abc8c99d3213247b3)
2. sdsl-lite -- 'git clone https://github.com/cosmo-team/sdsl-lite.git' (commit 9fa981958a9d2ddade12d083548f2b09939514fb)
3. stxxl -- 'git clone https://github.com/stxxl/stxxl' (commit 5b9663e6b769748f3b3d3a9a779b4b89e24d7a27)
4. tclap -- 'git clone https://github.com/eile/tclap' (commit f41dcb5ce3d063c9fe95623193bba693338f3edb)
5. Boost 1.54* -- 'wget http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.bz2'

* VARI fails to compile with later versions of Boost.  For the time being, it is necessary to download and compile Boost 1.54 and update BOOST_PATH in the Makefile to reflect the installed directory.  See Issue [#7](/../../issues/7).

They should be configured and built following their own instructions and set to install their files in a 3rd_party_inst subdirectory which is a sibling of 3rd_party_src.  The following sequence of commands should build the required parts.  Compilation errors may or may not affect the functionality of VARI, as VARI doesn't use all functionality of 3rd party sources.   Please email me if you run into trouble. I'm intermitently working on streamlining the process. -MDM May 17, 2017

**Note**: Change "/home/martin_muggli/git/test/cosmo" to wherever your cosmo working tree ends up.

    # Fetch software and setup directories
    git clone https://github.com/cosmo-team/cosmo/
    cd cosmo/
    git checkout VARI
    mkdir 3rd_party_src
    mkdir -p 3rd_party_inst/boost
    cd 3rd_party_src
    git clone https://github.com/refresh-bio/KMC
    git clone https://github.com/cosmo-team/sdsl-lite.git
    git clone https://github.com/stxxl/stxxl
    git clone https://github.com/eile/tclap
    wget http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.bz2
    tar -xjf boost_1_54_0.tar.bz2

    # Build the five dependencies
    cd boost_1_54_0
    ./bootstrap.sh --prefix=../../3rd_party_inst/boost
    ./b2 install
    cd ..
    
    cd sdsl-lite/
    /usr/bin/time sh install.sh /home/martin_muggli/git/test/cosmo/3rd_party_inst
    cd ..

    cd stxxl
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/martin_muggli/git/test/cosmo/3rd_party_inst -DBUILD_STATIC_LIBS=ON
    make
    make install
    cd ../..

    cd KMC
    make
    cd ..

    cd tclap/
    autoreconf -fvi
    ./configure --prefix=/home/martin_muggli/git/test/cosmo/3rd_party_inst
    make
    make install
    cd ..
    
    # Build VARI
    make


## Authors of VARI

Implemented by [Alex Bowe][abowe]. Original concept and prototype by [Kunihiko Sadakane][ksadakane].
Colored extension prototyped by Robert Raymond and substantially extended by Martin D. Muggli, with help from Alex Bowe.

These people also proved *incredibly* helpful: [Rayan Chikhi][rchikhi], [Simon Puglisi][spuglisi],
[Travis Gagie][tgagie], [Christina Boucher][cboucher], [Simon Gog][sgog], [Dominik Kempa][dkempa].


[dsk]: http://minia.genouest.org/dsk/
[minia]: http://minia.genouest.org/
[abyss]: https://github.com/bcgsc/abyss
[succ]: http://alexbowe.com/succinct-debruijn-graphs
[debby]: http://github.com/alexbowe/debby

[boost]: http://www.boost.org
[bgl]: http://www.boost.org/doc/libs/1_56_0/libs/graph/doc/
[sdsl-lite]: https://github.com/simongog/sdsl-lite
[networkx]: https://networkx.github.io/
[stxxl]: http://stxxl.sourceforge.net/
[python]: https://www.python.org/
[numpy]: http://www.numpy.org/
[tclap]: http://tclap.sourceforge.net/

[semver]: http://semver.org/
[nucleotides]: http://nucleotid.es/
[tci]: https://travis-ci.org

[abowe]: https://github.com/alexbowe
[cboucher]: http://christinaboucher.com/
[tgagie]: http://www.cs.helsinki.fi/u/gagie/
[ksadakane]: http://researchmap.jp/sada/
[spuglisi]: http://www.cs.helsinki.fi/u/puglisi/
[dkempa]: http://www.cs.helsinki.fi/u/dkempa/
[rchikhi]: https://github.com/rchikhi
[sgog]: https://github.com/simongog/


## Steps for running omGraph

1. Error correct optical maps using cOMet.

2. Build the de Bruijn graph using cosmo-build. Follow instructions from the cosmo-vari page. The output is a .dbg file.

https://github.com/cosmo-team/cosmo/tree/VARI

3. Use vari_rest_kmers to extract all restriction kmers from the graph

./vari_rest_kmers -o $RES_ENZ1 $KMER_SIZE.reads.dbg

4. Use vari_find_paths to find all simple paths between restriction nodes

./vari_find_paths <big_d> <factoring_loops> <dbg_file_path> <restriction_nodes_file_path>

ex
./vari_find_paths 100000 500 $KMER_SIZE.reads.dbg restriction_nodes

5. Use seeding_ext to find alignments of Rmaps in the DBG

./seeding_ext <simple_paths_file> <k_gram_size> <rmap_file_path>

ex 
./seeding_ext /ufrc/boucher/kingdgp/cosmo/cosmo/proper_edges/edges1en63_f500 3 ecoli.map


## Inputs of omGraph

1. De bruijn graph produced by cosmo. The extension is .dbg.

2. Rmap file in the following format. The fragments are in kbp.

rmap_0  1.259  5.397  6.57  2.696  22.841  17.412  1.836  3.224  1.51  10.925  9.356  11.331  8.78  7.887  50.162  28.192  8.706  5.601 
rmap_2  57.633  5.228  10.023  11.951  6.854  1.158  4.903  6.293  3.094  22.78  17.565  1.863  3.22  1.388  10.775  9.414  
rmap_3  5.935  6.568  2.709  22.6  17.239  1.823  3.198  1.508  10.965  9.313  11.361  8.8  7.966  49.833  27.797  8.571  5.577  
rmap_4  1.482  5.34  6.609  2.718  22.259  17.991  1.909  3.175  1.249  10.744  9.403  11.327  9.018  7.915  49.928  27.578  8.481 
rmap_5  5.688  6.491  2.705  22.663  17.396  1.931  3.261  1.412  10.957  9.437  11.33  8.89  7.767  50.386  27.165  13.687  


## Outputs of omGraph

1. The restriction nodes file contains the list of restriction nodes found in the graph along with their node ID in the graph

ex

54341 CGGACCGTTCACGCGCTTCAGCCTGTAAAAA

69190 CGGACCGATGAACTCTGGGTTCAGGCCAAAA

186848 CGGACCGTAATGCCTTACATTACCAGCCAAA

268644 CGGACCGGTCGTCGGGGAAATTCTCGAGAAA

294618 CGGACCGAAGAAACTGGCACCACCATCGAAA

308759 CGGACCGGCTAGTTACGCAGTGCGGCGGAAA

330034 CGGACCGCCAGGCGCAATTTGCACGATGAAA

336439 CGGACCGTCTCCATCTGGTCGTAGCCTGAAA


2. The simple paths file contains all simple paths between pairs of restriction nodes. 
The format is <source_node_id> <destination_node_id> <path_length_in_dbg>

ex

18954121 20846886 7730

18993618 16943406 12235

18993618 16998481 12235

18993618 16944255 12235

18993618 13568127 40501

19013547 3710101 16625

19211673 21148041 7133

19231575 10060163 4640

19292285 1160390 1951


3. The file rmap_alignements.txt contains the sequence of restriction nodes to which an Rmap aligns

rmap_0	640110 640123 22703333 8347638 15737788 15654806 640110 4985662 10816414 950755 

rmap_2	640110 640123 4985662 15956018 12597034 

rmap_4	12767312 640110 640123 10822804 15956018 12597034 

rmap_5	6221649 2671121 22333100 5260490 640123 640110 17178034 640110 4985662 10816414 

rmap_6	17736748 15180487 7818715 17515656 640123 640110 10923440 


4. The file rmap_alignements_lengths.txt contains the lengths of simple paths between an aligned pair of restriction nodes

ex

rmap_0	7936 50566 27997 8621 5656 0 35287 12150 12665 

rmap_2	7936 49741 27453 8762 

rmap_4	9403 7936 49441 27453 8762 

rmap_5	9252 11292 8692 7936 50566 26832 14093 52150 12150 

rmap_6	5656 8621 27997 50566 7936 8151 

(There is a simple path from 640110 to 640123 whose length is 7936)


