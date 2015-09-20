## Flexbar — flexible barcode and adapter removal

The program Flexbar preprocesses high-throughput sequencing data efficiently. It demultiplexes barcoded runs and removes adapter sequences. Moreover, trimming and filtering features are provided. Flexbar increases read mapping rates and improves genome as well as transcriptome assemblies. It supports next-generation sequencing data in fasta/q and csfasta/q format from Illumina, Roche 454, and the SOLiD platform. Contact [jtroehr](https://github.com/jtroehr) in case you need support with this SeqAn application.

### Reference

Matthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich: Flexbar — flexible barcode and adapter processing for next-generation sequencing platforms. Biology 2012, 1(3):895-905.

See article on [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24832523).

### Download

A version of Flexbar sources for [building](https://github.com/seqan/flexbar#building-from-source) can be obtained using git. The program is also provided as precompiled binary for Linux 64 and Mac OS X. Follow [installation](https://github.com/seqan/flexbar#installation-of-binaries) instructions after [download](https://github.com/seqan/flexbar/releases) of binaries. Additionally, Flexbar is available via package manager on Debian systems.

Older versions than Flexbar 2.4 can be found [here](https://sourceforge.net/projects/flexbar/files).

### Building from source

Make sure that `git` and `cmake` commands are available, as well as development and runtime files of the TBB library 4.0 or later (Intel Threading Building Blocks). Using a package manager is a simple way to install them. Furthermore, the [SeqAn](https://github.com/seqan/seqan) library 1.4.2 is required. Please [download](https://github.com/seqan/seqan/releases/download/seqan-v1.4.2/seqan-library-1.4.2.tar.bz2) it and decompress.

Get Flexbar source code:

        git clone https://github.com/seqan/flexbar.git Flexbar

Move SeqAn include folder to Flexbar:

        mv seqan-library-1.4.2/include Flexbar

Type these commands for building:

        cd Flexbar
        cmake .
        make

### Installation of binaries

To run Flexbar, the tbb library has to be available. Binary downloads contain the corresponding library file for runtime. Extract the downloaded archive and follow the platform specific instructions below.

#### Linux
Adjust lib search path to include the absolute path of the extracted Flexbar directory containing the lib file libtbb.so.2 for the current terminal session, or permanently in shell startup scripts:

        export LD_LIBRARY_PATH=/path/FlexbarDir:$LD_LIBRARY_PATH

#### Mac OS X
It applies the same as for Linux. Make the file libtbb.dylib available by setting the lib search path:

        export DYLD_LIBRARY_PATH=/path/FlexbarDir:$DYLD_LIBRARY_PATH

### Program features

* Demultiplexing of barcoded sequencing runs
* Detection and removal of adapter sequences
* Exact global alignment with free end-gaps
* Paired reads and separate barcode reads
* Wildcard N for barcodes and adapters
* Basic read filtering and trimming features
* Compressed input and output file support
* Extensive logging features, e.g. alignments
* Galaxy tool definition available in Tool Shed
* Multithreaded computation based on TBB library
* Sequence analysis based on SeqAn library

### Program usage

Flexbar needs at least one file with sequencing reads in fasta/q or csfasta/q format as input. Additionally, the target name, quality format of reads and further options can be specified. For read separation based on barcodes and for adapter removal, a file in fasta format with barcode or adapter sequences should be provided.

#### Synopsis

        flexbar -r reads [-t target] [-b barcodes] [-a adapters] [options]

Please refer to the help screen `flexbar -h` or [manual](https://github.com/seqan/flexbar/wiki).

#### Examples

        flexbar -r reads.fq -f i1.8 -t target -b brc.fa -a adap.fa

In this example, barcoded reads in illumina version 1.8 fastq format are demultiplexed by specifying a file with barcodes in fasta format. After read seperation based on barcodes, adapters given in fasta format are removed from the right side if they align at the read beginning or downstream. After removal the left side of reads is kept. Remaining reads are written to the file target.fastq in same format.

        flexbar -r reads.csfastq.gz -a adap.fa -ao 5 -ae LEFT -c

The second example, shows how to remove adapters in fasta format from left side of gzip compressed color-space (c) reads with quality scores (csfastq), if the overlap of adapter and read has at least length five. For left trim-end type the right side of reads is retained.

### Test data

To run Flexbar with the test datasets, make sure `flexbar` is reachable via the path variable and run `flexbar_validate.sh` within the test folder. Although default parameters of Flexbar are optimized to deliver good results in a large number of scenarios, the adjustment of parameters might improve results, e.g. `--adapter-min-overlap` and `--adapter-threshold`.

