#!/usr/bin/perl


######### Simulation ##########################################
#                                                             #
#       **************************************************    #  
#       **** Deisotope program for MS2		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################
package Spiders::Dta;

use strict;
use warnings;
use File::Basename;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
        _dta_file => undef,
    };
    bless $self, $class;
    return $self;
}
		
sub set_dta_file
{
	my ($self,$dtafile) = @_;
	return $self->{'_dta_file'} = $dtafile;	
}

sub get_dta_file
{
	my $self = shift;
	return $self->{'_dta_file'};
}

sub parse_dtafile_name
{
	my $self=shift;
	my $dtafile = $self->get_dta_file();
	$dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
	my ($scan, $specscan, $specscan2, $charge) = ($3, $3, $4, $5); 
	$self->{'_scan'} = $scan;
	$self->{'_charge'} = $charge;
}

sub get_scan_num
{
	my $self = shift;
	return $self->{'_scan'};
}

sub get_charge
{
	my $self = shift;
	return $self->{'_charge'};
}

sub process_dtafile
{
	my $self = shift;
	my $dtafile = $self->get_dta_file();

# define the mz hash for storing the mz and intensity
    my %mz_hash;
# open each dta file to load the data 
    open(DTAFILE,$dtafile) || die "can not open the dta file: $dtafile";
# get the precusor mz and charge 
    my $prec_mz_charge = <DTAFILE>;
    my ($prec_mz,$prec_charge) = split(/\s+/,$prec_mz_charge);
	my @mz_array;
	my @int_array;
    while(<DTAFILE>)
    {
# get the mz and intensity and save in the mz_hash
		my @data =split(/\s+/,$_);
		$mz_hash{$data[0]} = $data[1];
		push (@mz_array,$data[0]);
		push (@int_array,$data[1]);
	}
	$self->{'_precMZ'} = $prec_mz;
	$self->{'_mz_int_hash'} = %mz_hash;
	$self->{'_mz_array'}=\@mz_array;
	$self->{'_int_array'} = \@int_array;
}

sub get_prec_mz
{
	my $self = shift;
	return $self->{'_precMZ'};
}

sub get_mz_int_hash
{
	my $self = shift;
	return $self->{'_mz_int_hash'};
}

sub get_mz_array
{
	my $self = shift;
	return $self->{'_mz_array'};
}

sub get_int_array
{
	my $self = shift;
	return $self->{'_int_array'};
}

1;
