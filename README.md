## Flexbar — flexible barcode and adapter removal

The program Flexbar preprocesses high-throughput sequencing data efficiently. It demultiplexes barcoded runs and removes adapter sequences. Moreover, trimming and filtering features are provided. Flexbar increases read mapping rates and improves genome as well as transcriptome assemblies. It supports next-generation sequencing data in fasta and fastq format, e.g. from Roche 454 and the Illumina platform.

This is the Flexbar repository by Johannes Roehr. Flexbar is in the process of being adapted to SeqAn and incorporates features from the seqan flexcat repository. Refer to the [manual](https://github.com/seqan/flexbar/wiki) or contact [jtroehr](https://github.com/jtroehr) for support with this application.

### Reference

Matthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich: Flexbar — flexible barcode and adapter processing for next-generation sequencing platforms. Biology 2012, 1(3):895-905.

See article on [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24832523).

### Download

Flexbar source code as well as precompiled binaries for Linux and Mac OS X can be downloaded on the [release](https://github.com/seqan/flexbar/releases) page. Please follow instructions for building or installation of binaries below. Additionally, Flexbar is available via package manager on Debian systems. Versions before 2.4 can be found on the [old](https://sourceforge.net/projects/flexbar) page.

### Building from source

Make sure that `cmake` is available, as well as development and runtime files of the TBB library 4.0 or later (Intel Threading Building Blocks). Using a package manager is a simple way to install them. Furthermore, the SeqAn library is required:

* Get SeqAn library version 2.1.1 [here](https://github.com/seqan/seqan/releases)
* Download Flexbar source code [release](https://github.com/seqan/flexbar/releases)

Decompress both and move SeqAn include folder to Flexbar:

        mv seqan-library/include flexbar

Use these commands for building:

        cd flexbar
        cmake .
        make

For releases prior to 2.7 use SeqAn library 1.4.2 instead.

### Installation of binaries

To run Flexbar binaries after download, the corresponding TBB library has to be available. Downloads contain the library file for runtime. Follow the platform specific instructions below.

#### Linux
Adjust lib search path to include the absolute path of the Flexbar directory containing the lib file libtbb.so.2 for the current terminal session, or permanently in shell startup scripts:

        export LD_LIBRARY_PATH=/path/FlexbarDir:$LD_LIBRARY_PATH

#### Mac OS X
It applies the same as for Linux. Make the file libtbb.dylib available by setting the lib search path:

        export DYLD_LIBRARY_PATH=/path/FlexbarDir:$DYLD_LIBRARY_PATH

### Program usage

Flexbar needs at least one file with sequencing reads in fasta or fastq format as input. Additionally, the target name and further options can be specified. For read separation based on barcodes and for adapter removal, a file in fasta format with barcode or adapter sequences should be provided.

#### Synopsis

        flexbar -r reads [-b barcodes] [-a adapters] [options]

Please refer to the help screen `flexbar -h` or [manual](https://github.com/seqan/flexbar/wiki).

#### Examples

In this example, reads that are barcoded on left side are demultiplexed by specifying a file with barcodes in fasta format. After separation of reads, given adapters are removed from the right side if they do not align before read start. Subsequently, the left side of reads is kept if long enough. Remaining reads are written to the file `target.fastq` in same format as the input.

		flexbar -r reads.fq -t target -b brc.fa -be LEFT_TAIL -a adp.fa

### Test data

To run Flexbar with the test datasets, make sure `flexbar` is reachable via the path variable and run `flexbar_validate.sh` within the test folder. Although default parameters of Flexbar are optimized to deliver good results in many scenarios, the adjustment of parameters might improve results, e.g. `--adapter-min-overlap` and `--adapter-threshold`.

