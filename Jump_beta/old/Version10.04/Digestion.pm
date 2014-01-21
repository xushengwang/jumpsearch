package Spiders::Digestion;
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

require Exporter;
use Spiders::DatabaseModif;
use Spiders::Params;

$VERSION='0.02';
@ISA=qw(Exporter);
@EXPORT=qw(enzDigest readFasta writeFasta writeFreq);

########### version 0.02 ####################################################
# To reduce the memory usage, save the masshash, peptidehash and proteinhash#
# on the disc. Keep the filepoint on the memory                             #
#                                                                           #
#############################################################################

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


sub enzDigest{
	my ($self, $sqs, $paraHash) = @_;
	my $params = $self->get_parameter();
	my $enzyme = $paraHash->{'enzyme_info'};
	my $p = new Spiders::Params();
	my ($modif,$largest_modif) =$p->get_dynamic_modifications($params);	
	my $masshash;
	my $peptidehash;
	my $proteinhash;
	
##### $j is peptideID #######
	my $j = 0;

#### $i is proteinID #############
	for my $i (0..$#{$sqs}){
		my $id= $sqs->[$i]->{'id'};
		print "building index for protein: $id\n";
		my $desc = $sqs->[$i]->{'desc'};		
		
		my $peps = $self->digest($enzyme, $sqs->[$i]->{'seq'}, $paraHash->{'mis_cleavage'});
		foreach my $p (sort @{$peps}){
			
			next if (!defined ($p));
# changed by xusheng; remove the unknown amino acids			
			next if ($p =~ /X/);
			next if ($p =~ /U/);
			
#			$j++;
			next if(length($p)<4);
			
			$mw = $self->MW($p);
			next if ($mw < $paraHash->{'min_mass'} || $mw > $paraHash->{'max_mass'});
	################## partial tryptic ##############					
			$h->{'pepsGroup'} = $self->generatePeps($p) if ($paraHash->{'digestion'} eq 'partial');
	#################################################
	
			$mw = $mw * 1000;
			$mw = int($mw); 
#### remove these lines after add modified peptides 		
#			push @{$masshash->{$mw}},$j;
		
#			$peptidehash->{$j}->{'seq'}=$p;
#			$peptidehash->{$j}->{'proteinID'}=$i;
############################################################			
############# generate modified peptide #####################			
			if((scalar keys %$modif)>0)
			{
				my $modif_masshash;
				my @temp_peptide=();
				push (@temp_peptide,$p);
				$modif_masshash->{$mw} = \@temp_peptide;
				
				$self->generate_modif_peptide($modif_masshash);
				foreach my $mass ( keys(%$modif_masshash) ) 
				{
			###### if the mass of modified peptide larger than 6000, next
					next if ($mass > ($paraHash->{'max_mass'} * 1000));	
					
					foreach my $pep (@{$modif_masshash->{$mass}})
					{
						my $mass_mw = int($mass); 
						$j++;
						push @{$masshash->{$mass_mw}},$j;
						$peptidehash->{$j}->{'seq'}=$pep;
						$peptidehash->{$j}->{'proteinID'}=$i;
					}						
				}

			}
################# no dynamic modif ##########################			
			else
			{
				$j++;
				push @{$masshash->{$mw}},$j;
				$peptidehash->{$j}->{'seq'}=$p;
				$peptidehash->{$j}->{'proteinID'}=$i;				
			}
#########################################################			
		}
		$proteinhash->{$i}->{'id'} = $id;
		$proteinhash->{$i}->{'desc'} = $desc;
	}
	return ($masshash,$peptidehash,$proteinhash);
	
}

sub generate_modif_peptide
{
	my ($self,$pephash) = @_;
	my $cmdatabase = new Spiders::DatabaseModif();
	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");
	my $params = $self->get_parameter();
	my $combNUM = $params->{'max_modif_num'};
	my $p = new Spiders::Params();
	my ($modif,$largest_modif) =$p->get_dynamic_modifications($params);

	my %new_modif=();
	my $i=0;
	foreach $modif_aa (keys %$modif)
	{
		$new_modif{$modif_aa} = $modif_aa . $mod_symbol[$i];
		$i++;		
	}
	$cmdatabase->GenerateModifSeq($modif,$combNUM,$pephash,\%new_modif);
}

sub rmRedundancy{

	my ($self, $ps) = @_;
	my $h = {};
	foreach (@{$ps}){
		$h->{$_->{'id'}} = $_;
	}
	
	my $hh = [];
	foreach (sort keys %$h){
		push @{$hh},$h->{$_};
	}
	return $hh;
}


sub writeFreq{
	my ($self, $f, $freq) = @_;
	open OUT, '>'.$f or die;
	print OUT "Mass\tFrequency\n";
	foreach (keys %$freq){
		print OUT $_,"\t",$freq->{$_},"\n";
	}
	close OUT;
}



sub digest{

	my ($self, $enzyme, $seq, $NumOfMisCleavage) = @_;
	my @enzyme_array = split(/\s+/,$enzyme);
	my @sites = split(/\s*/,$enzyme_array[1]);
	my $uncleavage = $enzyme_array[2];
####################3 skip KP ##############
	for my $i (0..$#sites){
		$seq =~ s/($sites[$i])(?!$uncleavage)/$1 /g;
	}
	$seq =~ s/ \Z//;
	
	my @a = split(/ /,$seq);	
	
	if ($#a > -1){
	
		my @f = @a;
		
		for my $i (0..$#a){		
			my $ub = ($#a > $i+3) ? $i+3 : $#a;
			for my $j ($i+1..$ub){
				my @b = @a;
				push @f, join('',splice(@b,$i,$j-$i+1));
				 
			}
		}
		my $h = {};
		my $flag = 1;
		foreach my $frag (sort @f){			
			for my $i (0..$#sites){
				my $tmpFrag = $frag;
				my $site_uncleavage = $sites[$i] . $uncleavage;
				my $lc_site_uncleavage = lc($site_uncleavage);
				$tmpFrag =~ s/$site_uncleavage/$lc_site_uncleavage/;
				my $num = grep /$sites[$i]/,split(//,$tmpFrag);
				$flag *= ($num > $NumOfMisCleavage) ? 0 : 1;
			}
			
			if ($flag){
				$h ->{$frag}++;
			}
			$flag = 1;
		}
		
		return [sort keys %$h];
	}else{
		return [@a];
	}
}

sub MW{

	my ($self, $sq) = @_;
	my $parameter = $self->get_parameter();
	
	my $aaMW = {	
				'G' => 57.021464,	
				'D' => 115.02694,
				'A' => 71.037114,
				'Q' => 128.05858,
				'S' => 87.032029,
				'K' => 128.09496,
				'P' => 97.052764,
				'E' => 129.04259,
				'V' => 99.068414,
				'M' => 131.04048,
				'T' => 101.04768,
				'H' => 137.05891,
				'C' => 103.00919,
				'F' => 147.06841,
				'I' => 113.08406,
				'L' => 113.08406,
				'R' => 156.10111,
				'N' => 114.04293,	
				'Y' => 163.06333,
				'W' => 186.07931,				
				};
### add the static modification ################## 		
	foreach my $param (keys %$parameter)
	{
		if($param =~/add\_(\w+)\_/ and $parameter->{$param}>0)
		{
			$aaMW->{$1} += $parameter->{$param};
		}
	}	
			
	my $h = {};
	my @aas = split(//,$sq);
	foreach (@aas){
		$h->{$_}++;
	}
	my $mw = 0;
	foreach $aa (keys %$h){
		if(!defined($aaMW->{$aa}))
		{
#			print "\ryour sequences contain unkown amino acid: $aa!\n";
			$mw=0;
			last;
		}
		$mw += $h->{$aa} * $aaMW->{$aa};
	}
	### added the H20 and H+
	$mw += 19.017806;
	
	return $mw;
}

sub readFasta{

	my ($self, $f) = @_;
	
	my $parsedSqs = [];
	open F, $f or die;
	my @lines = <F>;
	close F;
	
	my @sqs = split("\n\>",join('',@lines));
	$sqs[0] =~ s/^\>//;
	
	
	for my $i (0..$#sqs){
		my $h = {};
		
		my @a = split("\n",$sqs[$i]);

		
		$h->{'id'} = substr($a[0], 0, index($a[0],' '));
		$h->{'desc'} = substr($a[0], index($a[0],' ')+1);
		shift @a;
		$h->{'seq'} = join('',@a);
		push @{$parsedSqs}, $h;
	}
	return $parsedSqs;
}

sub writeFasta{
	my ($self, $f, $sqs) = @_;
	open OUT, '>'.$f or die;
	for my $i (0..$#{$sqs}){
		print OUT '>'.$sqs->[$i]->{'id'},' ',$sqs->[$i]->{'desc'},"\n";
		print OUT $sqs->[$i]->{'seq'},"\n";
	}
	close OUT;
}

sub generatePeps{

	my ($self, $sq) = @_;
	
	my @a = split(//,$sq);
	
	my @b = split(//,reverse($sq));
	
	my @g;
	for my $i (0..$#a){
		my @t= @a;
		push @g, join('',splice(@t,0,$i+1));
	}
	for my $i (0..$#b){
		my @t= @b;
		push @g, join('',splice(@t,0,$i+1));
	}
	
	return [@g];

}

########## version 0.02 ######################




#############################################

1;

__END__