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


## Running omGraph

After the de Bruijn graph is built, the restriction nodes are identified using ./vari_rest_kmers. Then ./vari_find_paths finds the simple paths in the de Bruijn graph and aligns the optical maps to the paths. The file experiments/ecoli.sh gives the pipeline for the entire program. 
