#!/usr/bin/perl


######### Tag ##########################################
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

package Spiders::SummaryResults;

use strict;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new {
  my ($class,%args) = @_;
  my $self = {};
  $self->{PATH} = $args{-path};
  bless $self;
 
  return $self;
}

sub get_path
{
	my $self=shift;
	return $self->{PATH};
}

sub summary_tag
{
	my $self = shift;
	my $result_path = $self->get_path();
	my %Tag_num;
	my %freq;
	my @tag_files = glob("$result_path/*.tag");
	foreach (@tag_files)
	{
		open(TAGFILE,$_) || die "can not open the tag file!";
		my $number = scalar(<TAGFILE>);
		my $Tag_num{$_} = $number;
		if($number>10)
		{
			$number = 10;
		}
		$freq{$number}++;
	}
	return (\%Tag_num, \%freq);
}


sub Get_Tag_Peptide
{
	open(FILE,$ARGV[0]);
	my $string;
	print   "#   Rank/Sp      Id#     (M+H)+    deltCn   XCorr    Sp    Ions  Reference                           Peptide\tTag\trank_p\thyper_p\n";
	while(<FILE>)
	{
		chomp $_;
		if($_=~/^\s+1\..*\d+\s*\/\s*\d+.*\d*/)
		{
			my @array=split(/\s+/,$_);
			print $_;
			$string = $array[-1];
		}
		elsif($_=~/Tag\:\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/)
		{
			print "\t",$_,"\t";
			my $tag = $1;
			my $tag1 = $4; 

			$string =~ s/I/J/g;
			$string =~ s/L/J/g;
			
					$string =~ s/[^A-Z]//g;
			
			$tag=~ s/L/J/g;
			$tag1=~ s/L/J/g;
			$tag=~ s/[^A-Z]//g;
			$tag1=~ s/[^A-Z]//g;

			my $rev_tag = reverse ($tag);
			my $rev_tag1 = reverse ($tag1);
			
	#		print "aaaaaaaaa",$tag,"aaa\t";
	#		print "cccccc",$tag1,"ccc\t";
			
			if($string=~ /($tag)/)
			{
				print "match\t";
			}
			elsif($string=~ /($rev_tag)/)
			{
				print "match\t";
			}		
			else
			{
				print "unmatch\t";		
			}
			if($string=~ /($tag1)/)
			{
				print "match\t";
			}
			elsif($string=~ /($rev_tag1)/)
			{
				print "match\t";
			}		
			else
			{
				print "unmatch\t";		
			}
			
		}
		elsif($_=~/.*\.out/)
		{
			print "\n",$_;
		}
	}
}