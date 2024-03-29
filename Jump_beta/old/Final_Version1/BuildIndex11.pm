#!/usr/bin/perl

######### Simulation ##########################################
#                                                             #
#       **************************************************    #
#       **** Coverter MzXML to dta files	          ****    #
#       ****					                      ****    #
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #
#       ****all rights reserved.		              ****    #
#       ****xusheng.wang@stjude.org		              ****    #
#       ****					                      ****    #
#       ****					                      ****    #
#       **************************************************    #
###############################################################

package Spiders::BuildIndex;

use strict;
use warnings;
use File::Basename;
use Spiders::Digestion;
use Spiders::MakeDB;


use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
    };
    bless $self, $class;
	return $self;
}


sub set_parameter
{
	my ($self,$param)=@_;
	$self->{'_parameter'}=$param;
}


sub get_parameter
{
	my $self=shift;
	return $self->{'_parameter'};
}

sub get_parameter_hash
{
	my $self=shift;
	my $params = $self->get_parameter();
	my $fasta_file = $params->{'database_name'};
	my $enzyme_info = $params->{'enzyme_info'};
	my $digestion = $params->{'digestion'};	
	my $mis_cleavage = $params->{'max_mis_cleavage'};
	my $min_mass = $params->{'min_peptide_mass'};	
	my $max_mass = $params->{'max_peptide_mass'};	
			
	my $parameterHash = {'db_file'=>$fasta_file,
			'enzyme_info'=> $enzyme_info,
			'mis_cleavage'=> $mis_cleavage, # for example 2
		     'tryptic'=>$digestion,  # 'full|partial',
		     'min_mass'=>$min_mass,
	 	    'max_mass'=>$max_mass};
			
	return $parameterHash;
}

sub set_database
{
	my ($self,$database)=@_;
	$self->{'database'}=$database;
}


sub get_database
{
	my $self=shift;
	my $parameterHash = $self->get_parameter_hash();
	my $dta_path = $self->get_dta_path();
	my $data_file = $parameterHash->{'db_file'};	
	
	my $databasename = basename($data_file);
	$self->$self->{'database'} = "$dta_path/.database/$databasename";

	return $self->{'database'};
}


sub set_dta_path
{
	my ($self,$dtapath) = @_;

	if(!(-e $dtapath))
	{
		system(qq(mkdir $dtapath >/dev/null 2>&1));
	}
	$self->{_dta_path}=$dtapath;
}

sub get_dta_path
{
	my ($self) = shift;
	return $self->{_dta_path};
}


sub create_index {
	my ($self,$data_file) = @_;
	
	my $parameterHash = $self->get_parameter_hash();
	my $params = $self->get_parameter();
#	my $data_file = $parameterHash->{'db_file'};	

	
	#First create the index hashes 
#	my $parameterHash = {'mis_cleavage'=> 2,
#		     'digestion'=>'full',
#		     'min_mass'=>400,
#	 	    'max_mass'=>6000};
			

	my $makedb = new Spiders::MakeDB('db_file'=>$data_file);
	my $sqs = $makedb->make_db();
	
	my $digestion = new Spiders::Digestion();
	$digestion->set_parameter($params);	
	my ($masshash, $peptidehash,$proteinhash) = $digestion->enzDigest($sqs, $parameterHash);
	#Write the hashes results to an index file.
	return ($masshash, $peptidehash,$proteinhash);
}

sub create_index_file
{
	my ($self,$masshash,$peptidehash,$proteinhash,$databasename) = @_;
	my $dta_path = $self->get_dta_path();
	
#	my $databasename = basename($data_file);	
	if(! -e("$dta_path/.database/"))
	{
		system(qq(mkdir "$dta_path/.database/" >/dev/null 2>&1));
	}
	
	my $mass_index = $dta_path . "/.database/$databasename.mdx";
	my $peptide_index = $dta_path . "/.database/$databasename.pdx";
	my $protein_index = $dta_path . "/.database/$databasename.prdx";
	$self->CreateHashs($masshash,$peptidehash,$proteinhash,$mass_index,$peptide_index,$protein_index);
}

sub CreateHashs {
	my $self = shift;
	my $masshash = shift;	
	my $peptidehash = shift;
	my $proteinhash = shift;	
	my $masIndex = shift;
	my $peptideIndex = shift;
	my $proteinIndex = shift;
	
	#Check if the hash is not empty
	if(scalar(keys %{$masshash})<=0){ 
		die "Please verify that the file that you provided has a correct structure\n";
	}
	
	my @masses = sort {$a <=> $b} keys(%{$masshash}); #get the existing masses.	
	
	open(MOUTPUT,">>$masIndex") or die "Error opening file :$masIndex, $!\n";
	open(PEPOUTPUT,">>$peptideIndex") or die "Error opening file :$peptideIndex, $!\n";	
	foreach my $mas(0..$#masses){		

		print MOUTPUT pack("L L",$masses[$mas],$mas); #first write the mass
		#print MOUTPUT pack("N",$mas); #Write the mass id;
		
		#the petides will be organized by their mass weight first then by their sequence		
		my @pepIds =sort {$peptidehash->{$a}->{'seq'} cmp $peptidehash->{$b}->{'seq'}} @{$masshash->{$masses[$mas]}};
		
		foreach my $peptide (@pepIds){
			print PEPOUTPUT pack("L",$peptide); #PeptideID
			print PEPOUTPUT pack("L", $peptidehash->{$peptide}->{'proteinID'}); #ProteinID
			print PEPOUTPUT pack("L", $mas);#MassID
			print PEPOUTPUT pack("A255", $peptidehash->{$peptide}->{'seq'}); #Sequence
		}
	}
	close(MOUTPUT);
	close(PEPOUTPUT);
	
	open(PROOUTPUT,">>$proteinIndex") or die "Error opening file :$proteinIndex, $!\n";
	
	my @proteins = sort keys(%{$proteinhash});	
####### it must be sort, otherwise it is not in right order ############	
	foreach my $p (sort {$a<=>$b} @proteins){
		print PROOUTPUT pack("L",$p); #ProteinID
		print PROOUTPUT pack("A50",$proteinhash->{$p}->{"id"}); #Protein Name
		print PROOUTPUT pack("A200",$proteinhash->{$p}->{"desc"}); #Protein Description
	}

	close(PROOUTPUT);
}

sub readMassIndex{
	my ($self,$filename,$mass) = @_;
	
	my $fileSize = -s $filename;
	my $entrySize = 2*length(pack("L",0));
	my $nbElements = int($fileSize / $entrySize);		

	#Do a dicotomic search
	my $entry=$self->dicoSearch($filename,$nbElements,"L L",$entrySize,$mass,0);

	return $entry;		
}


sub getAssociatedPeptides{
	my ($self,$massId,$filename)= @_;
	
	#To search for a peptide first we do a dicotomic search on to find an peptide related to massID.
	#then we look into the elements next to the peptide to get the other peptides that have the same mass (if any)
	
	my $fileSize = -s $filename;
	my $entrySize = 3 * length(pack("L",0))+length(pack("A255","a"));
	my $nbElements = int($fileSize/$entrySize);
	
	#do dicotomic search and find a peptide 
	my ($entry,$index)= $self->dicoSearch($filename,$nbElements,"L L L A255",$entrySize,$massId,2);
	
	my @vals = unpack("L L L A255",$entry);
	my $peptidehash;
	$peptidehash->{$vals[0]}->{"seq"}= $vals[3];
	$peptidehash->{$vals[0]}->{"proteinID"}= $vals[1];
	
	#do a sequencial search 
	#first search on the left of the element
	my $curMassId = $massId;
	my $curIndex = $index-1;	
	my $input;
	open($input,"<$filename") or die "Error opening file $fileSize. $!\n";
	
	while($curMassId==$massId && $curIndex>0){
		my $val = $self->readEntry($input,$entrySize,$curIndex);	
		@vals = unpack("L L L A255",$val);
		
		if($vals[2]==$massId){
			$peptidehash->{$vals[0]}->{"seq"}= $vals[3];
			$peptidehash->{$vals[0]}->{"proteinID"}= $vals[1];
		}
		
		$curMassId = $vals[2];
		$curIndex--;
	}
	
	#Search on the right of the element
	$curMassId = $massId;
	$curIndex = $index+1;	
			
	while($curMassId==$massId && $curIndex<$nbElements){
		my $val = $self->readEntry($input,$entrySize,$curIndex);	
		@vals = unpack("L L L A255",$val);
		
		if($vals[2]==$massId){
			$peptidehash->{$vals[0]}->{"seq"}= $vals[3];
			$peptidehash->{$vals[0]}->{"proteinID"}= $vals[1];
		}
		
		$curMassId = $vals[2];
		$curIndex++;
	}
	
	close($input);
	
	return $peptidehash;
}

sub getProtein{
	my ($self,$proteinId,$filename) = @_;
	
	#Here we do a direct access;
	my $entrySize = length(pack("L",0))+length(pack("A50","x"))+length(pack("A200","x"));
	my $fileSize = -s $filename;
	
	#Here we suppose that the Ids are ordered sequencly , there is no missing IDs in the middle
	if($proteinId*$entrySize>$fileSize){
		die "the id you provided is out of range, please verify\n";
	}
	my ($input,$curIndex);
	open($input,"<$filename") or die "Error opening file $fileSize. $!\n";
	$curIndex = $proteinId;
	my $val = $self->readEntry($input,$entrySize,$curIndex);
	my @vals = unpack("L A50 A200",$val);	
	close($input);
	
	return ($vals[1],$vals[2]); #we return the ID and the description	
}

sub readEntry{
	my ($self,$file,$entrysize,$entryNum) = @_;
	my $entry;
	my $offset = $entrysize * ($entryNum);
	
	seek($file, $offset, 0) or return;
    read($file, $entry, $entrysize);

	return $entry;
}

sub dicoSearch{
	
	my ($self,$filename,$nbElements,$patern,$entrySize,$value,$valIndex) = @_;	
	my $lowBound =0;
	my $upBound = $nbElements-1;		
	my $input;
	open($input,"<$filename") or die "Error opening file $filename. $!\n";
	
	while(1){
		my $midPoint = int(($lowBound+$upBound)/2);
		my $midEntry = $self->readEntry($input,$entrySize,$midPoint);
		
		my @midvalue = unpack($patern,$midEntry);		
		
#		print $value,"\t",$midvalue[$valIndex],"aa\n";
		if($lowBound==$midPoint){return (-1,-1);}
		if($midvalue[$valIndex]==$value){						
			return ($midEntry,$midPoint);
		}
		else{ if($midvalue[$valIndex]>$value){
				$upBound = $midPoint;
			}
			else
			{
				$lowBound = $midPoint;
			}
		}				 			
	}
	close($input);
}

1;
