# Comet program #

This webpage is designed for testing Comet program. To run the Comet program, four files, including Comet source code, parameter file, human database, and a sample raw file, can be downloaded below:

Comet program can be downloaded [here](http://sourceforge.net/projects/comet-ms/files/comet_source.2014011.zip/download)

The parameter used for testing can be downloaded [here](https://code.google.com/p/jumpsearch/source/browse/wiki/comet.params?spec=svn34&r=34)

Human database can be downloaded [here](ftp://ftp.stjude.org/pub/software/jump/Human_20140709.fasta)

A testing raw file with mzXML can be downloaded [here](ftp://ftp.stjude.org/pub/software/jump/mm_15000.mzXML)

# Command used #

./comet.exe -Pcomet.params mm\_15000.mzXML

# Testing result #

All PSMs were sorted by E-value after searching. We detected 4805 and 6427 at 1% and 10% FDR, respectively.

# Linux System #
The program was run on Linux CentOS 6.5 operating system.