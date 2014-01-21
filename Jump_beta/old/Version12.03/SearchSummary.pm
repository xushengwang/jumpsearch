package Spiders::SearchSummary; 


######### SearchSummary ##########################################
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

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
        _result_dir => $arg{'result_dir'},
    };
    bless $self, $class;
    return $self;
}

sub get_result_dir
{
	my $self=shift;
	return $self->{'_result_dir'};
}

sub summary_spSearch
{
	my $self=shift;
	my $dir = $self->get_result_dir();
	
	my @spout = glob("$dir/*.spout");

	foreach my $spout (@spout)
	{
		open(SPOUT, $spout) || die "can not open the file:$spout\n";
		my $scannum = "";
		my $tag;
		my $tagEvalue;
		my $result;
		my %result;
		while(<SPOUT>)
		{
				if($_=~/scan number = (.*)/)
				{
						$scannum = $1;
				}
				if($_=~/Tag = (.*)\n/)
				{
						$tag = $1;
				}
				if($_=~/Tag E_value = (.*)\n/)
				{
						$tagEvalue = $1;
				}								
				if($_ =~/^\s+1\s+\d*\s+\d*.*\:\:/)
				{
					my @data = split(/\s+/,$_);
					$result{$data[4]}=$_;
					$result = $_;
				}
		}
		if(defined $result)
		{
			print $scannum,"\t",$tag,"\t",$tagEvalue,"\t",$result;
		}
=head		
		foreach $pep_score (sort {$result{$a}<=>$result{$b}} keys %result)
		{
				print $scannum,"\t",$result{$pep_score};;
				last;
		}
=cut
		
	}
}

sub summary_sequest
{
	my $self = shift;
	my $dir = $self->get_result_dir();
	
	my @out = glob("$dir/*.out");

	foreach my $out (@out)
	{
		open (IN, $out);
		$out =~ s/(\.)(\d+)(\.)(\d)(\.out)/$1$2$3$4$5/;
		my ($scan,$charge) = ($2,$4);

		seek(IN, 0,0);
		my ($mass, $database, $temp, $dCntemp, $dCntemp2, $dCntemp3, $dCntemp4) = grep(/\+\smass\s\=/ || /proteins\s\=/ || /^\s+1\.\s+\d+\s+\// || /^\s+2\.\s+\d+\s+\// || /^\s+3\.\s+\d+\s+\// || /^\s+4\.\s+\d+\s+\// || /^\s+5\.\s+\d+\s+\//, <IN>);

		next if (!defined($temp) || !defined($dCntemp));

			$temp =~ s/\s*\/\s*/\//g;
			$temp =~ s/(\s)\s+/$1/g;
			$temp =~ s/^\s+//; 
			$temp =~ s/\s+\Z//;

			$dCntemp =~ s/\s*\/\s*/\//g; 
			$dCntemp =~ s/(\s)\s+/$1/g; 
			$dCntemp =~ s/^\s+//; 
			$dCntemp =~ s/\s+\Z//;

			if (defined($dCntemp2)){
				$dCntemp2 =~ s/\s*\/\s*/\//g; $dCntemp2 =~ s/(\s)\s+/$1/g; $dCntemp2 =~ s/^\s+//; $dCntemp2 =~ s/\s+\Z//;
			if (defined($dCntemp3)){$dCntemp3 =~ s/\s*\/\s*/\//g; $dCntemp3 =~ s/(\s)\s+/$1/g; $dCntemp3 =~ s/^\s+//; $dCntemp3 =~ s/\s+\Z//;}
			if (defined($dCntemp4)){$dCntemp4 =~ s/\s*\/\s*/\//g; $dCntemp4 =~ s/(\s)\s+/$1/g; $dCntemp4 =~ s/^\s+//; $dCntemp4 =~ s/\s+\Z//;}
			}
			$mass =~ s/mass\s+\=\s+([\d\.]+)\s//;
			$mass = $1;
			$database =~ s/\s+//g; 
			my @temparray = split (/\,/, $database); 
			$database = pop(@temparray);

			if($database=~/hdr/)
			{
				$database=~s/\.hdr$//;	
			}
			

  
			my ($a, $rank, $id, $expmass, $c, $XCorr, $Sp, $e, $protein, $g, $h) = split(/\s/, $temp);
		 
			my ($z, $y,$i, $x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp);
			
			($z, $y,$i ,$x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp2) if ($dCn <= 0.01 && defined($dCntemp2));
			if (defined($dCntemp3)){($z, $y, $i, $x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp3) if ($dCn <= 0.01);}
			if (defined($dCntemp4)){($z, $y, $i, $x, $dCn, $w, $v, $u, $t, $s, $r) = split(/\s/, $dCntemp4) if ($dCn <= 0.01);}

			my $peptide = $g; 
			my $red = 0; 
			if(defined($h) ){
				$red = $g;
				$red =~ s/\+//;
				$peptide = $h;
			}
			my $ions = $e;
			$protein =~ s/\,.*//;
			print $out,"\t",$peptide,"\t",$XCorr,"\t",$Sp,"\t",$rank,"\t",$dCn,"\t",$protein,"\t",$mass,"\t",$expmass,"\n";
		}	
}

1;
