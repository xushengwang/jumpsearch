# JUMP v1.0 #

This is a software for database search.

## Introduction ##

JUMP is a new database search engine for peptide identification. JUMP incorporates sequence tag generation and tag scoring into database search.

## Usage ##

The basic usage of JUMP is quite simple:

JUMP.pl -param <JUMP parameter file> -rawfile < MS/MS data file>

> -param file              specifies a JUMP parameter file

> -rawfile  file           specifies an MS/MS data file

An MS/MS data file can be either .RAW or mzXML file that can be converted from .RAW file by various software, such as ReAdW or msconverter.

