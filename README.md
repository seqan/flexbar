## Flexbar — flexible barcode and adapter removal

The program Flexbar preprocesses high-throughput sequencing data efficiently. It demultiplexes barcoded runs and removes adapter sequences. Moreover, trimming and filtering features are provided. Flexbar increases read mapping rates and improves genome as well as transcriptome assemblies. It supports next-generation sequencing data in fasta and fastq format, e.g. from Illumina and the Roche 454 platform.

Refer to the [manual](https://github.com/seqan/flexbar/wiki) or contact [jtroehr](https://github.com/jtroehr) for support with this SeqAn application.

### Reference

Matthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich: Flexbar — flexible barcode and adapter processing for next-generation sequencing platforms. Biology 2012, 1(3):895-905.

See article on [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24832523).

### Download

Flexbar source code for building and precompiled binaries for Linux 64 and Mac OS X can be obtained on the release [download](https://github.com/seqan/flexbar/releases) page. Follow instructions for building or installation of binaries below. Versions older than 2.4 can be found [here](https://sourceforge.net/projects/flexbar). Additionally, Flexbar is available via package manager on Debian systems.

### Building from source

Make sure that `cmake` is available, as well as development and runtime files of the TBB library 4.0 or later (Intel Threading Building Blocks). Using a package manager is a simple way to install them. Furthermore, the SeqAn library is required.

* Download Flexbar source code release
* Get SeqAn library version 1.4.2 [here](https://github.com/seqan/seqan/releases/download/seqan-v1.4.2/seqan-library-1.4.2.tar.bz2)

Decompress both and move SeqAn include folder to Flexbar:

        mv seqan-library-1.4.2/include flexbar

Use these commands for building:

        cd flexbar
        cmake .
        make

### Installation of binaries

To run Flexbar binaries after download, the corresponding TBB library has to be available. Downloads contain the library file for runtime. Follow the platform specific instructions below.

#### Linux
Adjust lib search path to include the absolute path of the Flexbar directory containing the lib file libtbb.so.2 for the current terminal session, or permanently in shell startup scripts:

        export LD_LIBRARY_PATH=/path/FlexbarDir:$LD_LIBRARY_PATH

#### Mac OS X
It applies the same as for Linux. Make the file libtbb.dylib available by setting the lib search path:

        export DYLD_LIBRARY_PATH=/path/FlexbarDir:$DYLD_LIBRARY_PATH

### Program usage

Flexbar needs at least one file with sequencing reads in fasta or fastq format as input. Additionally, the target name, quality format of reads and further options can be specified. For read separation based on barcodes and for adapter removal, a file in fasta format with barcode or adapter sequences should be provided.

#### Synopsis

        flexbar -r reads [-t target] [-b barcodes] [-a adapters] [options]

Please refer to the help screen `flexbar -h` or [manual](https://github.com/seqan/flexbar/wiki).

#### Examples

        flexbar -r reads.fq -f i1.8 -t target -b brc.fa -a adap.fa

In this example, barcoded reads in illumina version 1.8 fastq format are demultiplexed by specifying a file with barcodes in fasta format. After separation of reads, adapters given in fasta format are removed from the right side if they do not align before read start. After removal the left side of reads is kept if long enough. Remaining reads are written to the file `target.fastq` in the same format.

		flexbar -r reads.fq.gz -a adap.fa -ao 5 -ae LEFT

The second example shows how to remove adapters in fasta format from left side of compressed reads with quality scores, if the overlap of adapter and read has at least length five. For left trim-end mode the right side of reads is retained.

### Test data

To run Flexbar with the test datasets, make sure `flexbar` is reachable via the path variable and run `flexbar_validate.sh` within the test folder. Although default parameters of Flexbar are optimized to deliver good results in a large number of scenarios, the adjustment of parameters might improve results, e.g. `--adapter-min-overlap` and `--adapter-threshold`.

