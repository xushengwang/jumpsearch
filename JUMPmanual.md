**JUMP v1.0**


**Introduction**


JUMP (Jumbo Mass
Spectrometry-based Proteomics Tool)
is a new database search engine for peptide identification. JUMP incorporates
sequence tag generation and tag scoring into database search.




**Usage**


The basic
usage of JUMP is quite simple:


JUMP.pl -p <parameter file>  <MS/MS data file>


-p parameter file             specifies a JUMP parameter file

MS/MS data file               specifies an MS/MS data file


An MS/MS data file can be
either .RAW or mzXML file that can be converted from .RAW file by various
software, such as ReAdW or msconverter.


**Parameters guide**
| **Name** | **Description** | **Type** | **Default**|
|:---------|:----------------|:---------|:-----------|
| database\_name | Specifies the FASTA or index protein database for search. If a FASTA  file was specified, the program will automatically build an indexed   database| string   | None       |
| enzyme\_info | Specifies the enzyme information used for creating a database  | string   | Tryptic KR P |
| digestion | Specifies the digestion methods. Three types of methods can be used,  including Full, Partial, and None | string   | Full       |
| max\_mis\_cleavage | A missed cleavage is defined as a site within the peptide matching  one of the cleavage rules (see enzyme\_info). This parameter is used to set  the maximum number of missed cleavages. | integer  | 2          |
| min\_peptide\_mass | Specifies the minimum mass value of peptides in the database | integer  | 400        |
| max\_peptide\_mass | Specifies the maximum mass value of peptides in the database | integer  | 6000       |
| inclusion\_decoy | This parameter is to set whether the database contains the decoy  sequences | boolean  | 1 (true) 0 (false) |
| isolation\_window | Specifies the isolation window used for peak selection, with unit of  m/z | float    | 1.0        |
| intrascanppm | ppm for peak matching within the same scan | integer  | 10         |
| interscanppm | ppm for peak matching between different scans | integer  | 5          |
| preproc\_msscanrange | range of ms scans used for decharging | integer  | 4          |
| ms2\_consolidation | Specifies number of peaks to keep within 100 m/z units   | integer  | 6          |
| ms2\_deisotope | Speficies whether to perform deisotope for MS2 peaks | boolean  | 1 (true) 0 (false) |
| prec\_mass\_window | Specifies to remove the peaks with the window with precursor mass +/-  value/2 (additional static 0.25 tolerance to either side) | float    | 3.0        |
| ppm      | Tolerance (ppm) used for MS2 deisotope  | integer  | 30         |
| M\_M1    | The intensity ratio between the peaks and following peaks.  | float    | 0.3        |
| Mn\_Mn1  | The intensity ratio between the peaks and following peaks. | Float    | 1          |
| charge12\_ppm | Specifies mass tolerance  | Integer  | 30         |
| ppi\_percentage |                 |          | 10         |
| tag\_tolerance | Mass tolerance used for tag inference | Float    | 0.02 (high) 0.5 (low) |
| low\_mass | The minimum mass allowed to define a residue | Integer  | 57         |
| high\_mass | The maximum mass allowed to define a residue | Integer  | 180        |
| tag\_select\_method | Two tag scoring methods can be used, including rank\_p and hyper\_p.  The rank\_p method calculates the score based on wilcox rank sum method; the  hyper\_p calculates the score based on hypergeometric method. | string   | rank\_p (high) hyper\_p (low) |
| ion\_series |  a, b, c, d, v, w, x, y, and z ions respectively. The values  entered for these paramters should either 0 or 1.  | String   |            |
| peptide\_tolerance | This parameter is used to find candidate peptides. The candidate peptides  is selected if the candidate mass is within this tolerance of the  experimental precursor mass. The units can be specified by the  ‘peptide\_tolerance\_units’ parameter. | integer  | 15         |
| peptide\_tolerance\_units | Two types of mass units: 1 for Da and 2 for parts per million (ppm). | integer  | 2          |
| frag\_mass\_tolerance | Specifies the fragment mass tolerance. The fragment tolerance should  be set at 0.01 for high depended on high-resolution instrument resolution; it  should be set at 0.5 for low-resolution instrument.  | float    | 0.01 (high) 0.5 (low) |
| tag\_search\_method | Two search methods can be specified: 1 for ; 2 for searching by the a  certain number of tags  | integer  | 1 or 2     |
| max\_number\_tags\_for\_search | This parameter is used to define the number of tags that are used for  search. The tag is sorted in terms of its p-value. The software search | integer  | 100        |
| tag\_coeff\_evalue | Specifies the weight of tag E value for calculating weighted E value. | float    | 1.0        |
| dynamic\_AA | This parameter is used to define dynamic amino acids.  | Float    |            |
| add\_AA\_name | This is static modification. If the value for any residues is setted,  it should always be treated as having a modification on their natural mass.  For example, the cysteine is carboxymethylated, this parameter would  be set to be: add\_C\_Cysteine = 57.02146. | float    |            |
| max\_modif\_num | This parameter is used to set the maximum number of dynamic  modifications that can be in the candidate peptide.  | integer  | 6          |
| cluster  | Specifies whether JUMP will be run on the cluster | boolean  | 1 (true) 0 (false) |
| Job\_Management\_System | Specifies the platform of the cluster. Three types of jobs can be  used: SGE, PBS and LSF. | String   | SGE        |
| number\_of\_selected\_result | Specifies the number of records being shown in the selected output section  in the output file.  | Integer  | 3          |
| number\_of\_detailed\_result | Specifies the number of records being shown in the “Selected  identified peptide” section in the output file. | integer  | 100        |
| simulation | Specifies whether to perform simulation.  If true, either one of  the following methods can used to generate a falsified data sets based on  input real data. | boolean  | 1 (true) 0 (false) |
| sim\_MS1 | If the ‘simulation’ parameter is set, the data set was generated by  adding precusor ion mass of each MS/MS spectrum by 100 ppm. | Integer  | 100        |
| sim\_MS2 | If the ‘simulation’ parameter is set, the data set was generated by  randomizing the MS/MS spectrum. Briefly, we generated uniform random number  from -5 to 5, and added this random number to each m/z of each peak. | integer  | 5          |


############################## Parameter file for JUMP
#################################


# JUMP version 10.02


# Date: 02/01/2013


# Author: Xusheng Wang


####################################################################################


############################## Database
#############################################


database\_name = /home/xwang4/JUMP\_database/human/human\_ft\_mc2\_c57.fasta.mdx


enzyme\_info = Tryptic KR P


digestion = full


max\_mis\_cleavage = 2


min\_peptide\_mass = 400.0000


max\_peptide\_mass = 6000.0000


inclusion\_decoy = 1


############################## Preprocessing ##########################################


precursor\_mass\_window = 0.19


intrascanppm = 10


interscanppm = 5


preproc\_msscanrange = 4


peptide\_tolerance = 15


peptide\_tolerance\_units = 2        #
1 = Th; 2 = ppm;


vary\_tolerance = 0


peptide\_tolerance\_frac = 0.8


ms2\_consolidation = 6


prec\_window = 3


ppm = 30


M\_M1 = 0.3


Mn\_Mn1 = 1


charge12\_ppm = 30


#############################
Tag##################################################


tag\_tolerance = 0.02


low\_mass = 57


high\_mass = 187


tag\_select\_method = rank\_p  # method can be input: hyper\_p or
rank\_p


######################Search
###########################################


ion\_series = 0 1 0 0 0 0 0 1 0  # a, b, c, d, v, w, x, y, and z
ions respectively


frag\_mass\_tolerance = 0.02


tag\_search\_method = 2       # 1. last when
found 2. search using the number of tag assigned by 'max\_number\_tag\_for\_search'


max\_number\_tags\_for\_search = 500  # If this number larger than the
total number of tag, the total number of tag will be used


tag\_coeff\_evalue = 4


##################### Dynamic Modification
######################################


# C:  57.02146 carbamidomethylation or 71.0371 acrylamide


# STY: 79.96633


# M: 15.99492


# GG: 114.04293


# SILAC K:4.02511, 6.02013, 8.01420


# SILAC R:6.02013, 10.00827


#STY = 0.000


#K = 0.000


dynamic\_M = 15.99492


max\_modif\_num = 6


#################### Static Modification
############################################


add\_Cterm\_peptide = 0.0000


add\_Nterm\_peptide = 0.0000


add\_G\_Glycine = 0.0000


add\_A\_Alanine = 0.0000


add\_S\_Serine = 0.0000


add\_P\_Proline = 0.0000


add\_V\_Valine = 0.0000


add\_T\_Threonine = 0.0000


add\_C\_Cysteine = 57.02146


add\_L\_Leucine = 0.0000


add\_I\_Isoleucine = 0.0000


add\_X\_LorI = 0.0000


add\_N\_Asparagine = 0.0000


add\_O\_Ornithine = 0.0000


add\_B\_avg\_NandD = 0.0000


add\_D\_Aspartic\_Acid = 0.0000


add\_Q\_Glutamine = 0.0000


add\_K\_Lysine = 0.0000


add\_Z\_avg\_QandE = 0.0000


add\_E\_Glutamic\_Acid = 0.0000


add\_M\_Methionine = 0.0000


add\_H\_Histidine = 0.0000


add\_F\_Phenylalanine = 0.0000


add\_R\_Arginine = 0.0000


add\_Y\_Tyrosine = 0.0000


add\_W\_Tryptophan = 0.0000


add\_J\_user\_amino\_acid = 0.0000


add\_U\_user\_amino\_acid = 0.0000


#################### Cluster system
####################################################


cluster = 1


Job\_Management\_System = SGE


#################### Simulation ######################################################


simulation = 0   # 1 = yes; 0 = no


sim\_MS1 = 0    # unit is ppm


sim\_MS2 = 0     # unit is Da


###################### Output for results
###############################################


number\_of\_selected\_result = 3


number\_of\_detailed\_result = 100


#####################################################################################


JUMP manual


**System Requirement**


JUMP must meet these minimum requirements.


Cluster settings:

Hardware:
·         Run on a cluster system

·         SGE or LSF job management system

·         32 GB memory on each node

·         Run on a single server

·         2 GHz processor

Software:

·         WINE install if run .raw file


·         Perl modules:

**Features**


JUMP includes the following features:


1.    Support for input
formats: .raw or .mzXML


2.    Support for FASTA
databases or indexed databases

continue to working on it

