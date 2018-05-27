## Flexbar – flexible barcode and adapter removal

The program Flexbar preprocesses high-throughput sequencing data efficiently. It demultiplexes barcoded runs and removes adapter sequences. Moreover, trimming and filtering features are provided. Flexbar increases read mapping rates and improves genome as well as transcriptome assemblies. Unique molecular identifiers can be extracted in a flexible way. The program supports sequencing data in fasta and fastq format, e.g. from the Illumina platform.

Refer to the [manual](https://github.com/seqan/flexbar/wiki) or contact [Johannes Roehr](https://github.com/jtroehr) for support with this application.

![Flexbar logo](https://github.com/seqan/flexbar/wiki/images/flexbar-logo.png)


### References

Johannes T. Roehr, Christoph Dieterich, Knut Reinert:  
Flexbar 3.0 – SIMD and multicore parallelization. Bioinformatics 2017.

See article on [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28541403)

Matthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich:  
Flexbar – flexible barcode and adapter processing for next-generation sequencing platforms. Biology 2012.

See article on [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/24832523)


### Download

Flexbar source code as well as binaries for Linux and Mac OS can be downloaded on the [release](https://github.com/seqan/flexbar/releases) page. Please follow instructions for building or setup of binaries below. Additionally, Flexbar is available via package manager on Debian systems. Versions before 2.4 can be found on the [old](https://sourceforge.net/projects/flexbar) page.


### Building from source

Make sure that `cmake` is available, as well as development and runtime files of the TBB library 4.0 or later (Intel Threading Building Blocks). For example on Debian systems, install the packages `libtbb-dev` and `libtbb2`. Furthermore, the SeqAn library and a compiler that supports C++14 is required:

* Get SeqAn library version 2.2.0 [here](https://github.com/seqan/seqan/releases/download/seqan-v2.2.0/seqan-library-2.2.0.tar.xz)
* Download Flexbar 3.2 source code [release](https://github.com/seqan/flexbar/releases)

Decompress both files:

	tar xzf flexbar-3.2.0.tar.gz
	tar xJf seqan-library-2.2.0.tar.xz

Move SeqAn include folder to Flexbar:

	mv seqan-library-2.2.0/include flexbar-3.2.0

Use these commands for building:

	cd flexbar-3.2.0
	cmake .
	make

Flexbar version 2.7 requires SeqAn 2.1.1 instead. Releases prior to 2.7 use the SeqAn 1.4.2 library.


### Binaries

For execution of provided Flexbar binaries, the corresponding TBB library has to be available. Downloads contain the library file for runtime. Follow the platform specific instructions below.

#### Linux
Adjust lib search path to include the absolute path of the Flexbar directory containing the lib file libtbb.so.2 for the current terminal session, or permanently in shell startup scripts:

	export LD_LIBRARY_PATH=/path/FlexbarDir:$LD_LIBRARY_PATH

#### Mac OS
It applies the same as for Linux. Make the file libtbb.dylib available by setting the lib search path:

	export DYLD_LIBRARY_PATH=/path/FlexbarDir:$DYLD_LIBRARY_PATH


### Program usage

Flexbar needs at least one file with sequencing reads in fasta or fastq format as input. Additionally, the target name and further options can be specified. For read separation based on barcodes and for adapter removal, a file in fasta format with barcode or adapter sequences should be provided.

	flexbar -r reads [-b barcodes] [-a adapters] [options]

Refer to the help screen `flexbar -h` or [manual](https://github.com/seqan/flexbar/wiki) for more information. Although default parameters of Flexbar are optimized to deliver good results in many scenarios, the adjustment of parameters might improve results, e.g. `--adapter-min-overlap`. To run tests, make sure `flexbar` is reachable via the path variable and run `flexbar_test.sh` within the test folder.

#### Examples

In this example, reads in fastq format are trimmed based on their quality scores in Illumina version 1.8 format. The TAIL method trims the right end of reads until a quality score equal or higher than the threshold is reached, default 20. The option `--qtrim-threshold` can be used to change the score threshold. Remaining reads are written to the file `target.fastq` in same format as the input.

	flexbar -r reads.fq -t target -q TAIL -qf i1.8

Reads that are barcoded on the left end are demultiplexed by specifying a file with barcodes in fasta format. Reads that can be assigned to barcodes are written to separate files using file names that are based on the names of barcodes in the fasta file.

	flexbar -r reads.fq -b barcodes.fa -be LTAIL

To remove adapter sequences from single-end reads, specify a file with adapters in fasta format. These are removed from the right side of reads per default, if they do not align before the read start. The left side of reads is kept if long enough. The overlap of an adapter and read must have at least length 3 with at most 20% errors in default settings.

	flexbar -r reads.fq -a adapters.fa -ao 3 -at 0.2

For paired-end libraries, specify both files with paired reads and a fasta file with adapters for removal. Given adapters are trimmed in right mode per default. Additionally, it is recommended to activate the overlap detection in case of normal paired reads. This increases the sensitivity by removing also very short parts of adapters if an overlap between paired reads is detected.

	flexbar -r r1.fq -p r2.fq -a adapters.fa -ap ON

For further examples visit the [manual](https://github.com/seqan/flexbar/wiki) page.

