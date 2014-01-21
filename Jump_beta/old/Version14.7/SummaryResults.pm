#!/usr/bin/perl


######### Tag #################################################
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
use List::Util qw(first);
use File::Basename;
#use Statistics::Descriptive;

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

######################
# Example 1: Extract tag rank information 
# my $folder = $ARGV[0];
# use Spiders::SummaryResults;
# my $ss = new Spiders::SummaryResults();
# my $mainhash=$ss->parse_spOut_files_v5($folder);
# foreach my $scan (keys %$mainhash)
#{
#	print $scan,"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"TagRank"},"\n";
#}
#####################################


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
	my $Tag_num;
	my %freq;
	my $freq_target_decoy;
	my @tag_files = glob("$result_path/*.tag");
	foreach my $file (@tag_files)
	{
		my $spout = $file;
		$spout =~ s/tag/spout/;
		if(-e "$result_path/$spout")
		{
			my $mainhash = $self->read_spout($spout);
			my $target_decoy = $self->check_decoy_target($mainhash);
			open(TAGFILE,"$file") || die "can not open the tag file!";
			my @lines = <TAGFILE>;
			my $number = scalar @lines;

			my $Tag_num->{$file} = $number;
			if($number>15)
			{
				$number = 10;
			}
			$freq{$number}++;
			if($target_decoy eq "Target")
			{
				$freq_target_decoy->{$number}->{"Target"}++;
			}
			else
			{
				$freq_target_decoy->{$number}->{"Decoy"}++;			
			}
		}
	}
	return ($Tag_num, \%freq,$freq_target_decoy);
}


sub get_spOut_files
{
	my $self = shift;
	my $result_path = $self->get_path();
	my @spfiles = glob("$result_path/*.tag");
	return \@spfiles;
}


sub checkParameters{
	my ($self,$i,@line) = @_;
	
	if($line[1]!~/\d+/ ){
		die "Error parsing file at line $i\n The MH field should be a numeric value\n";
	}
	
	if($line[2]!~/\d+/){
		die "Error parsing file at line $i\n The lSideMass field should be a numeric value\n";
	}
	
	if($line[3]!~/\d+/){
		die "Error parsing file  at line $i\n The rSideMass field should be a numeric value\n";
	}
	
	if($line[4]!~/\d+/){
		die "Error parsing file  at line $i\n The peptide E value  field should be a numeric value\n";
	}
}

sub check_decoy_target
{
	my ($self,$mainhash) = @_;
	
	my $target_decoy="";
	
	#Get files associated with each threshold keys	
	foreach my $file(keys %{$mainhash}){	
		if($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"reference"} =~ /Decoy/)
		{
			$target_decoy = "Decoy";
		}
		else
		{
			$target_decoy = "Target";		
		}
	}
	
	return $target_decoy;	
}

sub parse_spOut_files
{
	my ($self)= @_;
	my $folder = $self->get_path();	
	my $mainhash;
	my $maxHitCondiered = 10;
	
	if(!-e $folder){ die "Error opening folder $folder\n";}
	
	my @files = glob("$folder/*.spout");
	
	foreach my $spout (@files){

#		print "\rProcessing file $spout\n";		
		open(INPUT,$spout) or die "Could onpen $spout. Error:!$\n";
		my @content = <INPUT>;		
		$spout = basename($spout);
		$spout =~ s/\.spout//;		
		
		if(scalar(@content)>0)
		{ 
			($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}) = grep(/Precursor mass\s*=/,@content);
			($mainhash->{$spout}->{"outhash"}->{"PrecursorMatchedPep"}) = grep(/Precursor matched peptides\s*=/,@content);			
			($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}) = grep(/Tag number\s*=/,@content);
			
			chomp($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"});
			$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];
			$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/\s+/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];
			
			
			chomp($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"});
			$mainhash->{$spout}->{"outhash"}->{"TotalTagNum"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}))[1];
			
			my $startIndex = first {$content[$_]=~/-+\n/} 0..$#content;

			foreach my $i ($startIndex+1..$#content)
			{					
				last if ($content[$i]=~/^\s*$/);
				last if ($content[$i]=~/All identified peptide/);
				chomp($content[$i]); #Remove return line
				$content[$i]=~ s/^\s+//; #Remove leading spaces
				my @vals = split /\s+/,$content[$i];
				if(scalar(@vals)==1)
				{
					next;
				}
				if(scalar(@vals)==10)
				{
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"TagRank"} =0;
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"reference"}= $vals[8];	
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Weight_E_value"} = $vals[6];					
				}
			#	print $spout,"\t",scalar(@vals),"\n";
				elsif(scalar(@vals)==11){
				#	die "Error parsing file $spout at l:qine: $i\n Some values are missing\n";
				
				
					if($vals[0]!~/\d+/){
						die "Error parsing file $spout at line $i\n The Oder field should be a numeric value\n";
					}
				#	last if ($vals[0]>$maxHitCondiered);
				#	next if ($vals[9]=~/\@/);
					
					$self->checkParameters($i,@vals);
					print $spout,$vals[8],"\n"  if ($vals[8]<0.1);
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"MH"}= $vals[1];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"lsideMass"}= $vals[2];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"rsideMass"}= $vals[3];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide_E_value"}= $vals[4];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag"} = $vals[5];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"TagRank"} = $vals[6];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag_num"} = $vals[7];				
					
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Weight_E_value"} = $vals[8];

								
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"reference"}= $vals[9];
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide"}= $vals[10];

					my $intpep = $vals[10];
					$intpep =~ s/[A-Z\-]\.([A-Z\@\#\*]+)\.[A-Z\-]/$1/;
					
					$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"intpep"}= $intpep;
									
			#		$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"modification"}= $vals[7];
				}
			}
		}	 
		close($spout);		
	}
	return $mainhash;
}

sub Output_summary_table
{
	my ($self,$mainhash,$outfile) = @_;
	open(SUM,">$outfile") || die "can not open the summary file: $outfile";
	foreach (keys %$mainhash)
	{
        print SUM $_,"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"peptide"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"reference"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"peptide_E_value"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"Tag"},"\n";
	}
	close(SUM);
}

sub Output_tag_summary_table
{
	my ($self,$mainhash,$outfile) = @_;
	open(SUM,">$outfile") || die "can not open the summary file: $outfile";
	foreach (keys %$mainhash)
	{
        print SUM $_,"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"peptide"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"reference"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"peptide_E_value"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"Tag"},"\n";
	}
	close(SUM);
}



sub score_distribution
{
	my ($self,$mainhash) = @_;
	my $score;
	foreach my $file(keys %{$mainhash})
	{
		next if(!defined($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"Weight_E_value"}));
		my $pep_score = sprintf("%.1f",($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"Weight_E_value"}+0.2));
		
		if($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"reference"} =~ /Decoy/)
		{
	
			$score->{"Decoy"}->{$pep_score}++;
		}
		else
		{
			$score->{"Target"}->{$pep_score}++;			
		}
	}
	return $score;
}

sub FDR_calculation
{
	my ($self,$mainhash) = @_;
	
	my $Evals;
	my $PeptideEvals;
	my $scores_FDR_hash;
	
	#Get files associated with each threshold keys	
	foreach my $file(keys %{$mainhash}){	
		push @{$Evals->{int($mainhash->{$file}->{"outhash"}->{"Tag_E_value"}+0.5)}},$file;
		my $Weight_E_value = sprintf("%.1f", $mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"Weight_E_value"});
#		push @{$PeptideEvals->{int($mainhash->{$file}->{"outhash"}->{"Order"}->{1}->{"Weight_E_value"}+0.5)}},$file;	
		push @{$PeptideEvals->{$Weight_E_value}},$file;
	}
	
#	$scores_FDR_hash =$self->CalcTagEval($mainhash,$scores_FDR_hash,$Evals);
	$scores_FDR_hash= $self->CalcPeptideEvalue($mainhash,$scores_FDR_hash,$PeptideEvals);
	
	return $scores_FDR_hash;
}

sub CalcTagEval{
	my ($self,$mainhash,$scores_FDR_hash,$Evals)=@_;
	
	my $DecoyRefs_key;
	my $NonDecoryRefs_key;
	my @EvalKeys = sort {$a <=> $b} keys %{$Evals};
	my ($nbDecoyRef,$nbNonDecoyRef)=(0,0); #used to calculate FDR
	my ($nbSup,$nbTotal)=(0,0); #Used to calculate Coverage
############
		my $test_number = scalar keys %{$mainhash};

###############	
	foreach my $key (reverse sort {$a<=>$b} @EvalKeys){ #Go through the files associated with each Tag Evalue 
#################### changed 			
#		$nbTotal = scalar(@{$Evals->{$key}});		
#############################		
		foreach my $spout(@{$Evals->{$key}}){
													
			if($mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"reference"}=~/##Decoy__/){
				$nbDecoyRef++;
			}
			else
			{
				$nbNonDecoyRef++;
			}
			
			if($mainhash->{$spout}->{"outhash"}->{"Tag_E_value"}>$key){
				$nbSup++;
			}					 			 			 
		}
		
		 #Save how many ##Decoy and non ##Decoy refrences we have for each key, to use it when calculating
		 #the sensitivity and specificity
		 
		 $DecoyRefs_key->{$key} = $nbDecoyRef;
		 $NonDecoryRefs_key->{$key} = $nbNonDecoyRef;			 
		 			
		 if($nbNonDecoyRef>0){ 
			$scores_FDR_hash->{"Tag_E_value"}->{"$key"}->{"FDR"} = $nbDecoyRef / $nbNonDecoyRef;										
		 }
		 else #in case for one key all the files have only ##Decoy top refs.
		 {
			 $scores_FDR_hash->{"Tag_E_value"}->{"$key"}->{"FDR"} =-$nbDecoyRef;
		 }

############## changed #################################
		 $scores_FDR_hash->{"Tag_E_value"}->{"$key"}->{"nbDecoy"} = $nbDecoyRef;
		 $scores_FDR_hash->{"Tag_E_value"}->{"$key"}->{"nbTaget"} = $nbNonDecoyRef;
		 
		 $scores_FDR_hash->{"Tag_E_value"}->{"$key"}->{"Coverage"}= ($nbDecoyRef+$nbNonDecoyRef)/$test_number;
#		 $scores_FDR_hash->{"Tag_E_value"}->{"$key"}->{"Coverage"}= $nbSup/$nbTotal;
######################################
	}

	#Calculate the sensitivity and specificity for each threshold
	my @thresholds = sort {$a <=> $b} keys %{$DecoyRefs_key};
	foreach my $key (0...$#thresholds)
	{

		my $UpkeysDecoy =0;
		my $DownkeysDecoy =0;
		my $UpKeyNonDecoy=0;
		my $DownNonDecoy=0;
	
		for(my $i=$key;$i<=$#thresholds;$i++){
			$UpkeysDecoy += $DecoyRefs_key->{$i};
			$UpKeyNonDecoy += $NonDecoryRefs_key->{$i};			
#			$UpkeysDecoy= $DecoyRefs->{$thresholds[$i]};
#			$UpKeyNonDecoy= $NonDecoryRefs->{$thresholds[$i]};	
		}
		
		#Get number of ##Decoy and non #Decoy for keys < threshold
		for(my $i=$key-1;$i>=0;$i--){
			$DownkeysDecoy += $DecoyRefs_key->{$i};
			$DownNonDecoy += $NonDecoryRefs_key->{$i};
#			$DownkeysDecoy= $DecoyRefs->{$thresholds[$i]};
#			$DownNonDecoy= $NonDecoryRefs->{$thresholds[$i]};	
		}
		

		if($UpKeyNonDecoy+$DownNonDecoy!=0){		
			$scores_FDR_hash->{"Tag_E_value"}->{"$thresholds[$key]"}->{"sensitivity"} = $UpKeyNonDecoy/ ($UpKeyNonDecoy+$DownNonDecoy);
		}
		else
		{
			$scores_FDR_hash->{"Tag_E_value"}->{"$thresholds[$key]"}->{"sensitivity"}=0;
		}
		if($DownkeysDecoy+$UpkeysDecoy!=0){
			$scores_FDR_hash->{"Tag_E_value"}->{"$thresholds[$key]"}->{"specificity"} = $DownkeysDecoy / ($DownkeysDecoy+$UpkeysDecoy);
		}
		else
		{
			$scores_FDR_hash->{"Tag_E_value"}->{"$thresholds[$key]"}->{"specificity"} =0;
		}		
	}	
	return $scores_FDR_hash;
}

sub CalcPeptideEvalue{
	
	my ($self,$mainhash,$scores_FDR_hash,$PeptideEvals)=@_;
	
	my $DecoyRefs_key;
	my $NonDecoryRefs_key;
	my @EvalKeys = sort {$a <=> $b} keys %{$PeptideEvals};
	my ($nbDecoyRef,$nbNonDecoyRef)=(0,0); #used to calculate FDR (cumulative)
	my ($nbSup,$nbTotal)=(0,0); #Used to calculate Coverage (cumulative)

	my $test_number = scalar keys %{$mainhash};
	foreach my $key (reverse sort {$a<=>$b} @EvalKeys)
	{ #Go through the files associated with each Tag Evalue 

		my $nbDecoyRef_key = 0;
		my $nbNonDecoyRef_key = 0;
		
		$nbTotal += scalar(@{$PeptideEvals->{$key}});		

		foreach my $spout(@{$PeptideEvals->{$key}}){
			next if(!defined($mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"reference"}));
			if($mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"reference"}=~/##Decoy__/){
				$nbDecoyRef++;
				$nbDecoyRef_key++;
			}
			else
			{
				$nbNonDecoyRef++;
				$nbNonDecoyRef_key++;
			}
			
			if($mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"Weight_E_value"}>$key)
			{
				$nbSup++;
			}					 			 			 
		}
	
		 #Save how many ##Decoy and non ##Decoy refrences we have for each key, to use it when calculating
		 #the sensitivity and specificity

		 $DecoyRefs_key->{$key} = $nbDecoyRef_key;
		 $NonDecoryRefs_key->{$key} = $nbNonDecoyRef_key;			 
		 			
		 if($nbNonDecoyRef>0){ 
			$scores_FDR_hash->{"Weight_E_value"}->{"$key"}->{"FDR"} = $nbDecoyRef / $nbNonDecoyRef;										
		 }
		 else #in case for one key all the files have only ##Decoy top refs.
		 {
			 $scores_FDR_hash->{"Weight_E_value"}->{"$key"}->{"FDR"} =-$nbDecoyRef;
		 }			 
############## changed #################################
		
		 $scores_FDR_hash->{"Weight_E_value"}->{"$key"}->{"nbDecoy"} = $nbDecoyRef;
		 $scores_FDR_hash->{"Weight_E_value"}->{"$key"}->{"nbTaget"} = $nbNonDecoyRef;		 
		 $scores_FDR_hash->{"Weight_E_value"}->{"$key"}->{"Coverage"}= ($nbDecoyRef+$nbNonDecoyRef)/$test_number;
#		 $scores_FDR_hash->{"peptide_E_value"}->{"$key"}->{"Coverage"}= $nbSup/$nbTotal;
#		print $scores_FDR_hash->{"peptide_E_value"}->{"$key"}->{"nbTaget"},"\n";
######################################		 
		 
	}

	###### Total Decoy number ##########
	my $total_decoy = 0;
	foreach my $decoy_score (keys %{$DecoyRefs_key})
	{
		$total_decoy += $DecoyRefs_key->{$decoy_score};
	}
	print $total_decoy,"\n";
	my $total_target = 0;
	foreach my $target_score (keys %{$NonDecoryRefs_key})
	{
		$total_target += $NonDecoryRefs_key->{$target_score};
	}
	print $total_target,"\n";	

	#Calculate the sensitivity and specificity for each threshold
	my @thresholds = sort {$a <=> $b} keys %{$DecoyRefs_key};
	
	foreach my $key (0...$#thresholds){
		print $thresholds[$key],"\t";
		my $UpkeysDecoy =0;
		my $DownkeysDecoy =0;
		my $UpKeyNonDecoy=0;
		my $DownNonDecoy=0;

		#Get number of ##Decoy and non #Decoy for keys >= threshold
		for(my $i=$key;$i<=$#thresholds;$i++){
			$UpkeysDecoy += $DecoyRefs_key->{$thresholds[$i]};

			$UpKeyNonDecoy += $NonDecoryRefs_key->{$thresholds[$i]};			
#			$UpkeysDecoy= $DecoyRefs->{$thresholds[$i]};
#			$UpKeyNonDecoy= $NonDecoryRefs->{$thresholds[$i]};	
		}

		#Get number of ##Decoy and non #Decoy for keys < threshold
		for(my $i=$key-1;$i>=0;$i--){

			$DownkeysDecoy += $DecoyRefs_key->{$thresholds[$i]};
			$DownNonDecoy += $NonDecoryRefs_key->{$thresholds[$i]};
#			$DownkeysDecoy= $DecoyRefs->{$thresholds[$i]};
#			$DownNonDecoy= $NonDecoryRefs->{$thresholds[$i]};	
		}
		
		print $UpKeyNonDecoy,"\t",$DownNonDecoy,"\t",$DownkeysDecoy,"\t",$UpkeysDecoy,"\n";	
	
		if($UpKeyNonDecoy+$DownNonDecoy!=0){		
			$scores_FDR_hash->{"Weight_E_value"}->{"$thresholds[$key]"}->{"sensitivity"} = $UpKeyNonDecoy/ $total_target;
		}
		else
		{
			$scores_FDR_hash->{"Weight_E_value"}->{"$thresholds[$key]"}->{"sensitivity"}=0;
		}
		if($DownkeysDecoy+$UpkeysDecoy!=0)
		{
				$scores_FDR_hash->{"Weight_E_value"}->{"$thresholds[$key]"}->{"specificity"} = $UpkeysDecoy / $total_decoy;
		}
		else{
				$scores_FDR_hash->{"Weight_E_value"}->{"$thresholds[$key]"}->{"specificity"} =0;
			}
		
	}
	
	return $scores_FDR_hash;
}



sub XcorrFDR_calculation{
	my ($self,$runhash) = @_;
	
#	my $path = $self->get_path();
#	my $runhash = $self->parse_outfile($path);
	
	my $Xcorr;
	my $scores_FDR_hash;
	
	#Get files associated with each threshold keys	
	foreach my $file (keys %{$runhash}){
		my $xcorr_value = sprintf("%.2f", $runhash->{$file}->{"XCorr"});
		push @{$Xcorr->{$xcorr_value}},$file;
	}
	
	$scores_FDR_hash =$self->CalcXcorr($runhash,$scores_FDR_hash,$Xcorr);
	
	return $scores_FDR_hash;
}

sub CalcXcorr{
	
	my ($self,$runhash,$scores_FDR_hash,$Xcorr)=@_;
	
	my $DecoyRefs_key;
	my $NonDecoryRefs_key;
	my @Xcorr = sort {$a <=> $b} keys %{$Xcorr};
	my ($nbDecoyRef,$nbNonDecoyRef)=(0,0); #used to calculate FDR (cumulative)
	my ($nbSup,$nbTotal)=(0,0); #Used to calculate Coverage (cumulative)

	my $test_number = scalar keys %{$runhash};
	foreach my $key (reverse sort {$a<=>$b} @Xcorr)
	{ #Go through the files associated with each Tag Evalue 
		my $nbDecoyRef_key = 0;
		my $nbNonDecoyRef_key = 0;
		
		$nbTotal += scalar(@{$Xcorr->{$key}});		

		foreach my $out (@{$Xcorr->{$key}}){
			next if(!defined($runhash->{$out}->{"protein"}));		
			if($runhash->{$out}->{"protein"}=~/Decoy__/){
				$nbDecoyRef++;
				$nbDecoyRef_key++;
			}
			else
			{
				$nbNonDecoyRef++;
				$nbNonDecoyRef_key++;
			}
			
			if($runhash->{$out}->{"XCorr"}>$key){
				$nbSup++;
			}					 			 			 
		}
	
		 #Save how many ##Decoy and non ##Decoy refrences we have for each key, to use it when calculating
		 #the sensitivity and specificity

		 $DecoyRefs_key->{$key} = $nbDecoyRef_key;

		 $NonDecoryRefs_key->{$key} = $nbNonDecoyRef_key;			 
		
		 if($nbNonDecoyRef>0){ 
			$scores_FDR_hash->{"XCorr"}->{"$key"}->{"FDR"} = $nbDecoyRef / $nbNonDecoyRef;										
		 }
		 else #in case for one key all the files have only ##Decoy top refs.
		 {
			 $scores_FDR_hash->{"XCorr"}->{"$key"}->{"FDR"} =-$nbDecoyRef;
		 }			 
############## changed #################################
		
		 $scores_FDR_hash->{"XCorr"}->{"$key"}->{"nbDecoy"} = $nbDecoyRef;
		 $scores_FDR_hash->{"XCorr"}->{"$key"}->{"nbTaget"} = $nbNonDecoyRef;		 
		 $scores_FDR_hash->{"XCorr"}->{"$key"}->{"Coverage"}= ($nbDecoyRef+$nbNonDecoyRef)/$test_number;
######################################		 
		 
	}

	###### Total Decoy number ##########
	my $total_decoy = 0;
	foreach my $decoy_score (keys %{$DecoyRefs_key})
	{
		$total_decoy += $DecoyRefs_key->{$decoy_score};
	}
	print $total_decoy,"\n";
	my $total_target = 0;
	foreach my $target_score (keys %{$NonDecoryRefs_key})
	{
		$total_target += $NonDecoryRefs_key->{$target_score};
	}
	print $total_target,"\n";	

	#Calculate the sensitivity and specificity for each threshold
	my @thresholds = sort {$a <=> $b} keys %{$DecoyRefs_key};
	
	foreach my $key (0...$#thresholds){
		print $thresholds[$key],"\t";
		my $UpkeysDecoy =0;
		my $DownkeysDecoy =0;
		my $UpKeyNonDecoy=0;
		my $DownNonDecoy=0;

		#Get number of ##Decoy and non #Decoy for keys >= threshold
		for(my $i=$key;$i<=$#thresholds;$i++){
			$UpkeysDecoy += $DecoyRefs_key->{$thresholds[$i]};

			$UpKeyNonDecoy += $NonDecoryRefs_key->{$thresholds[$i]};			
#			$UpkeysDecoy= $DecoyRefs->{$thresholds[$i]};
#			$UpKeyNonDecoy= $NonDecoryRefs->{$thresholds[$i]};	
		}

		#Get number of ##Decoy and non #Decoy for keys < threshold
		for(my $i=$key-1;$i>=0;$i--){

			$DownkeysDecoy += $DecoyRefs_key->{$thresholds[$i]};
			$DownNonDecoy += $NonDecoryRefs_key->{$thresholds[$i]};
#			$DownkeysDecoy= $DecoyRefs->{$thresholds[$i]};
#			$DownNonDecoy= $NonDecoryRefs->{$thresholds[$i]};	
		}
		
		print $UpKeyNonDecoy,"\t",$DownNonDecoy,"\t",$DownkeysDecoy,"\t",$UpkeysDecoy,"\n";

		if($UpKeyNonDecoy+$DownNonDecoy!=0){		
			$scores_FDR_hash->{"XCorr"}->{"$thresholds[$key]"}->{"sensitivity"} = $UpKeyNonDecoy/ $total_target;
		}
		else
		{
			$scores_FDR_hash->{"XCorr"}->{"$thresholds[$key]"}->{"sensitivity"}=0;
		}
		if($DownkeysDecoy+$UpkeysDecoy!=0)
		{
			$scores_FDR_hash->{"XCorr"}->{"$thresholds[$key]"}->{"specificity"} = $UpkeysDecoy / $total_decoy;
		}
		else
		{
			$scores_FDR_hash->{"XCorr"}->{"$thresholds[$key]"}->{"specificity"} =0;
		}		
	}
	
	return $scores_FDR_hash;
}



sub summary_second_peptide
{
	my ($self,$folder)= @_;
	my $mainhash;
	
	if(!-e $folder){ die "Error opening folder $folder\n";}
	my $mainhash;
	my @files = glob("$folder/*.spout");
	foreach my $spoutfile (@files)
	{
		$spoutfile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).spout)\Z/$1/;
		my $spout = $3;
		my $order = $4;
		open(INPUT,$spoutfile) or die "Could onpen $spout. Error:!$\n";
		my @content = <INPUT>;		
		close($spoutfile);				
		if(scalar(@content)>0)
		{ #if the file is not empty
	#		($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"},$mainhash->{$spout}->{"outhash"}->{"Tag"},$mainhash->{$spout}->{"outhash"}->{"Tag_E_value"},$mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"}) = grep(/Precursor mass\s*=/ || /Tag\s*=/ || /Tag E_value\s*=/ || /Tag Side Mass\s+=/,@content);
			($mainhash->{$spout}->{$order}->{"outhash"}->{"PrecursorMass"}) = grep(/Precursor mass\s*=/,@content);
			
			chomp($mainhash->{$spout}->{$order}->{"outhash"}->{"PrecursorMass"}); $mainhash->{$spout}->{$order}->{"outhash"}->{"PrecursorMass"} = (split(/=/,$mainhash->{$spout}->{$order}->{"outhash"}->{"PrecursorMass"}))[1];
	#		chomp($mainhash->{$spout}->{"outhash"}->{"Tag"}); $mainhash->{$spout}->{"outhash"}->{"Tag"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"Tag"}))[1];
	#		chomp($mainhash->{$spout}->{"outhash"}->{"Tag_E_value"}); $mainhash->{$spout}->{"outhash"}->{"Tag_E_value"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"Tag_E_value"}))[1];
	#		chomp($mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"}); $mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"}))[1];
			
			my $startIndex = first {$content[$_]=~/-+\n/} 0..$#content;
					
			foreach my $i ($startIndex+1..$#content)
			{
				last if ($content[$i]=~/^\s*$/);
				last if ($content[$i]=~/All identified peptide/);
				chomp($content[$i]); #Remove return line
				$content[$i]=~ s/^\s+//; #Remove leading spaces
				my @vals = split /\s+/,$content[$i];
				if(scalar(@vals)==1)
				{
					next;
				}
				
				if(scalar(@vals)!=10){
					die "Error parsing file $spout at line: $i\n Some values are missing\n";
				}
				
				if($vals[0]!~/\d+/){
					die "Error parsing file $spout at line $i\n The Oder field should be a numeric value\n";
				}
				
				$self->checkParameters($i,@vals);
			

				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"MH"}= $vals[1];
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"lsideMass"}= $vals[2];
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"rsideMass"}= $vals[3];
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide_E_value"}= $vals[4];
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag"} = $vals[5];
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag_E_value"} = $vals[6];
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"Tag_Side_Mass"} = $vals[7];
							
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"reference"}= $vals[8];
				$mainhash->{$spout}->{$order}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"peptide"}= $vals[9];
		#		$mainhash->{$spout}->{"outhash"}->{"Order"}->{"$vals[0]"}->{"modification"}= $vals[7];
			}
		}	

	}
	return $mainhash;
	
}

sub TagDistribution
{
	my ($self,$mainhash,$file) = @_;
	open(FILE,">$file") || die "can not open the file: $file\n";
	open(SUMMARY,">sum_$file");
	my %target_tag_length;
	my %decoy_tag_length;
	foreach my $spout (keys %$mainhash)
	{
		if($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"} eq "N/A")
		{
			print FILE $mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"reference"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"},"\t","0","\n";
			if($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"reference"} =~ /Decoy/)
			{
				$decoy_tag_length{0}++;
			}
			else
			{
				$target_tag_length{0}++;
			}
		}
		elsif(defined($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"}))
		{
			print FILE $mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"reference"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"},"\t",length($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"}),"\n";
			if($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"reference"} =~ /Decoy/)
			{
				$decoy_tag_length{length($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"})}++;
			}
			else
			{
				$target_tag_length{length($mainhash->{$spout}->{"outhash"}->{"Order"}->{"1"}->{"Tag"})}++;
			}			
		}
	}
	foreach my $key (sort {$a<=>$b} keys %target_tag_length)
	{
		print $key,"\t",$target_tag_length{$key},"\t",$decoy_tag_length{$key},"\n";
	}
}



1;
