#!/usr/bin/perl -I /home/xwang4/scripts

package Spiders::Postprocessing;

######### Deisotope ##########################################
#                                                             #
#       **************************************************    #  
#       **** Postprocessing for search		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 2012 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################


use Spiders::RankHits;


require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw();
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.01;


sub get_mis_cleavage_number
{
	my ($self,$peptide) = @_;
	
	my $mis = ($peptide =~ s/([KR][A-OQ-Z])/$1/g) || 0;
	return $mis;
} 

sub get_mod_site_number
{
	my ($self,$peptide) = @_;
	my $mods = ($testpeptide =~ s/([\@\#\*\^\~\$]+)/$1/g) || 0;
	return $mods;
}

sub parse_spOut_files_v5{
	my ($self,$folder)= @_;
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
		{ #if the file is not empty
	#		($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"},$mainhash->{$spout}->{"outhash"}->{"Tag"},$mainhash->{$spout}->{"outhash"}->{"Tag_E_value"},$mainhash->{$spout}->{"outhash"}->{"Tag_Side_Mass"}) = grep(/Precursor mass\s*=/ || /Tag\s*=/ || /Tag E_value\s*=/ || /Tag Side Mass\s+=/,@content);
			($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}) = grep(/Precursor mass\s*=/,@content);
			($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}) = grep(/Tag number\s*=/,@content);
			
			chomp($mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}); $mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];

			$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"} = (split(/\s+/,$mainhash->{$spout}->{"outhash"}->{"PrecursorMass"}))[1];
			chomp($mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}); $mainhash->{$spout}->{"outhash"}->{"TotalTagNum"} = (split(/=/,$mainhash->{$spout}->{"outhash"}->{"TotalTagNum"}))[1];
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








sub calculate_mass_shift
{
	my ($mass, $expmass, $charge) = ($$outhash{'MH'}, $$outhash{'expMH'}, $$outhash{'charge'});
	my $realdiff = $mass-$expmass; 
    my $ppm = ($diff/$mass/$charge)*(10**6);
	
    if ($$outhash{'XCorr'} >= $XCorr)
	{
		next if ($$outhash{'protein'} =~ /Decoy__/);
		%{$goodhash{$outfile}} = %$outhash;
	}
	my ($sum, $sum2, $n) = (0,0,scalar(keys %goodhash));
	while (my ($outfile, $hash) = each %goodhash)
	{
		my $ppm = $$hash{'ppm'};
		$sum += $ppm; 
		$sum2 += $ppm**2;
	}
	my $premean= $sum/$n;
	my $prestdev = $utils->calc_stdev($sum, $sum2, $n);	
    if ($n < $bf){
		while (my ($outfile, $outhash) = each %$runhash){
			$$outhash{status} = 1;
		}
        print "Skipping accurate mass filtering because of not enough good scans.\n";
        return;
    }	
	my @vgarray;
	for (my $cycle=1; $cycle<=2; $cycle++)
	{
		my @array;
		($sum, $sum2, $n) = (0,0,0);
		for my $outfile (sort {$goodhash{$a}{'scan'}<=>$goodhash{$b}{'scan'}} keys %goodhash)
		{
			my $ppm = $goodhash{$outfile}{'ppm'};
			if (($ppm < $premean - 2*$prestdev) || ($ppm > $premean + 2*$prestdev))
			{
				$$runhash{$outfile}{'status'} = -1;
				delete ($goodhash{$outfile});
			}
			else 
			{
				push (@array, $outfile);
				$sum += $ppm;
				$sum2 += $ppm**2;
				$n++;
			}
		}
		$premean = $sum/$n;
		$prestdev = $utils->calc_stdev($sum, $sum2, $n);
		@vgarray = @array;
	}
	my $ppmshift = $sum/$n;
	
	return $ppmshift;	
}

sub calculate
{

}

sub calculate_FDR
{

	my $score = $ss->score_distribution($mainhash);

	my %merge_score;
	foreach (keys %{$score->{"Decoy"}})
	{
			$merge_score{$_}=1;
	}

	foreach (keys %{$score->{"Target"}})
	{
			$merge_score{$_}=1;
	}

	open(FDRRESULT,">FDR_Escore.txt");
	print FDRRESULT "E_score\tFDR\n";
	my $sum_target = 0;
	my $sum_decoy = 0;
	my $flag = 0;
	my $FDR_1_score = -1;
	my $prev_score = 0;
	foreach (reverse sort {$a<=>$b} keys %merge_score)
	{
			$score->{"Target"}->{$_}=0 if(!defined($score->{"Target"}->{$_}));
			$score->{"Decoy"}->{$_}=0 if(!defined($score->{"Decoy"}->{$_}));
			next if ($_==0);

			if($_<=100)
			{
					next if ($score->{"Target"}->{$_} == 0);
					my $FDR = $score->{"Decoy"}->{$_}/$score->{"Target"}->{$_};
					print FDRRESULT $_,"\t",$FDR,"\n";
			}
			$sum_target = $score->{"Target"}->{$_}+$sum_target;
			$sum_decoy = $score->{"Decoy"}->{$_}+$sum_decoy;
			if(($flag == 0) and ($sum_decoy*2/($sum_target+$sum_decoy)>0.01))
			{
					$FDR_1_score = $_;
					$flag = 1;
			}
			print RESULT2 $_,"\t",$sum_target,"\t",$sum_decoy,"\n";
			$prev_score = $_;
	}
}

my $tryptic = $utils->istryptic($peptide);

sub istryptic{  #returns 3 on fully tryptic, 2 on C-terminal partially tryptic, 1 on N-terminal partially tryptic, and 0 on non-tryptic
        shift @_;
        my ($Sequence) = @_;

        my $Nterm = 0;
        my $Cterm = 0;
        if ($Sequence =~ m/^[RK]\.[A-OQ-Z].*/ || $Sequence =~ m/^[\-]\..*/) {$Nterm = 1}; # Fixed bug DMD 3/20/2007
        if (($Sequence =~ m/.\.[A-Z\#\@\*]*[KR][\#\@\*]*\.[A-OQ-Z]\Z/) || ($Sequence =~ m/.\.[A-Z\#\@\*]*\.\-\Z/)) {$Cterm = 1};
        if ($Nterm && $Cterm) {
                return 3;
        } elsif ($Nterm || $Cterm){
                return 1 if $Nterm;
                return 2 if $Cterm;
        } else {
                return 0;
        }
}


sub change_score
{

}

