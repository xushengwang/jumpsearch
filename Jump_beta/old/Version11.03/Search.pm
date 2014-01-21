package Spiders::Search;

######### Simulation ##########################################
#                                                             #
#       **************************************************    #
#       **** Search database	                      ****    #
#       ****					                      ****    #
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #
#       ****all rights reserved.		              ****    #
#       ****xusheng.wang@stjude.org		              ****    #
#       ****					                      ****    #
#       ****					                      ****    #
#       **************************************************    #
###############################################################

use Spiders::Params;
use Spiders::MassUtils;
use Spiders::Hypergeometric;
use Spiders::BuildIndex;
use Spiders::MathUtils;


require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(findTag readTags exportMatches);
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.02;

######### version 1.02 ############
# To speed the search algorithm, it is changed to use the index file to search the database 
######### version 1.04 ############
# change the search module to allow search for low resolution data 
# 1. added the fragmented ion mass tolerance
######## Version 1.05 ##############
#add the amino acid to both sides of peptide  
######################################


sub new{

	my ($class, %arg) = @_;
	my $self = {};
	bless $self, $class;
	return $self
}

sub set_H_value
{
	my ($self,$c_value)=@_;
	$self->{_H_value}=$c_value;	
}

sub get_H_value
{
	my $self=shift;
	if(!defined($self->{_H_value}))
	{
		$self->{_H_value}=1.007276466812;
	}
	return $self->{_H_value};
}

sub readTags{

	my ($self, $f) = @_;
	
	my $tagHash = {};
	open F, $f or die;	
	while(<F>){
		if ($. > 1){
			chomp;
			# 6 columns: (1) Scan number, (2) precursor mass; (3) Tag sequences; (4) Side mass; (5) Rank p value; (6) Hyper p value 
			my @cols = split(/\t/,$_);
			my $h = {'scanNum'=>$cols[0],
					'precMass'=>$cols[1],
					'tagSeq'=>$cols[2],
					'sideMass'=>$cols[3],
					'rankP'=>$cols[4],
					'hyperP'=>$cols[5]
					};
			$tagHash->{$cols[0]} = $h;
		}
		
	}
	close F;

	return $tagHash;	
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

sub set_tag
{
	my ($self,$tag)=@_;
	$self->{'tag'}=$tag;
}

sub get_tag
{
	my ($self)=@_;
	return $self->{'tag'};
}

sub set_tag_number
{
	my ($self,$tag_num)=@_;
	$self->{'tag_num'}=$tag_num;
}

sub get_tag_number
{
	my ($self)=@_;
	return $self->{'tag_num'};
}

sub set_precMass
{
	my ($self,$precmass) = @_;
	$self->{'precMass'} = $precmass;
}

sub get_precMass
{
	my ($self) = @_;
	return $self->{'precMass'};
}

sub set_precCharge
{
	my ($self,$prec_charge) = @_;
	$self->{'precCharge'} = $prec_charge;
}

sub get_precCharge
{
	my ($self) = @_;
	return $self->{'precCharge'};
}

sub get_scanNum
{
	my $self=shift;
	my $tag=$self->get_tag();
	$self->{'scanNum'} = $tag->{'scanNum'};	
	return $self->{'scanNum'};
}

sub get_tagSeq
{
	my $self=shift;
	my $tag=$self->get_tag();
	$self->{'tagSeq'} = $tag->{'tagSeq'};	
	return $self->{'tagSeq'};
}

sub get_sideMass
{
	my $self=shift;
	my $tag=$self->get_tag();
	$self->{'sideMass'} = $tag->{'sideMass'};
	return $self->{'sideMass'};
}

sub set_database
{
	my ($self,$database) = @_;
	$database=~s/\.mdx//;
	$self->{'database'} = $database;
}

sub get_database
{
	my ($self) = @_;
	return $self->{'database'};
}

sub set_masshash
{
	my ($self,$masshash) = @_;
	$self->{'masshash'} = $masshash;
}

sub get_masshash
{
	my ($self) = @_;
	return $self->{'masshash'};
}

sub set_peptidehash
{
	my ($self,$peptidehash) = @_;
	$self->{'peptidehash'} = $peptidehash;
}

sub get_peptidehash
{
	my ($self) = @_;
	return $self->{'peptidehash'};
}

sub set_proteinhash
{
	my ($self,$proteinhash) = @_;
	$self->{'proteinhash'} = $proteinhash;
}

sub get_proteinhash
{
	my ($self) = @_;
	return $self->{'proteinhash'};
}


sub get_tagPvalue
{
	my $self=shift;
	my $tag=$self->get_tag();
	$self->{'rankP'} = $tag->{'rankP'};
	return $self->{'rankP'};
}

sub set_mass_tolerance
{
	my ($self,$mass_tolerance)=@_;
	$self->{'mass_tolerance'} = $mass_tolerance;
}

sub get_mass_tolerance
{
	my $self=shift;
	return $self->{'mass_tolerance'};
}

sub set_mass_tolerance_units
{
	my ($self,$mass_tolerance_units)=@_;
	$self->{'mass_tolerance_units'} = $mass_tolerance_units;
}

sub get_mass_tolerance_units
{
	my $self=shift;
	return $self->{'mass_tolerance_units'};
}



sub set_frag_tolerance
{
	my ($self,$mass_tolerance)=@_;
	$self->{'frag_tolerance'} = $mass_tolerance;
}

sub get_frag_tolerance
{
	my $self=shift;
	return $self->{'frag_tolerance'};
}

sub SearchMass{

	my ($self) = @_;
	my %MatchedMass;
	
#	my %proteinid;
#	my %precursor;	
#	my %tagGroups;
#	my %threotical_MH;
	
	my $index = new Spiders::BuildIndex();
	my $database = $self->get_database();
	
	my $mass_index_file = $database . ".mdx";
	my $peptide_index_file = $database . ".pdx";	
	my $protein_index_file = $database . ".prdx";	

	my $error = $self->get_mass_tolerance();
	my $error_unit = $self->get_mass_tolerance_units();
	
	my $precmass = $self->get_precMass();
#	my $H = $self->get_H_value();	
	
#	my $max_modif = 0;
######### version 1.1 using index to search the value #####################	
# get precursor mass 
# if there is no modification 
####### ----------------------o-------------------------o------------
####### ---------------------start---------------------end------------
#############################=============selected=======############
##### Version 1.12 both ppm and th can be used ################
	if($error_unit==2)
	{
		$error = $error*$precmass/1000000;
	}

	my $prec_start = int(($precmass-$error)*1000);
	my $prec_end = int(($precmass+$error)*1000);

#	$self->search_core(\%proteinid,\%tagGroups,\%threotical_MH,$prec_start,$prec_end);

	for(my $prec_mass = $prec_start; $prec_mass <= $prec_end; $prec_mass++)
	{
		my $massId = $index->readMassIndex($mass_index_file,$prec_mass);

		if($massId!=-1)
		{
			my $peptidehash =$index->getAssociatedPeptides($massId,$peptide_index_file);

			foreach $PeptideID (keys %{$peptidehash})
			{	
						
				my $mod;

				my %mod_pos=();
				
				my $nomod_pep = $peptidehash->{$PeptideID}->{'seq'};
				my @seq_array = split(/\./,$nomod_pep);
				$nomod_pep = $seq_array[1];

				while($nomod_pep =~ /[^a-zA-Z]+/g)
				{
					$mod_pos{$-[0]}=1;
					$nomod_pep =~ s/[^a-zA-Z]+//;
				}				
				if((scalar keys %mod_pos)==0)
				{
					my $length=length($peptidehash->{$PeptideID}->{'seq'});
					#$nomod_pep=$peptidehash->{$PeptideID}->{'seq'};
					$mod = ":" x $length;
				}
				else
				{
############## if the peptide has modification, get the modification info ##################	
# Fix a bug: K.QIVWKYCGR.M@			
#					my $nomod_pep = $peptidehash->{$PeptideID}->{'seq'};
#					my @seq_array = split(/\./,$nomod_pep);
#					$nomod_pep = $seq_array[1];			
					my @nomodseq = split(//,$nomod_pep);
				
					my $length=length($nomod_pep);
					for(my $i=0;$i<$length;$i++)
					{
						if($mod_pos{$i+1})
						{
							$mod .= ":";
							$mod .= $nomodseq[$i];
						}
						else
						{
							$mod .= ":";
						}
					}
				}

#				print $peptidehash->{$PeptideID}->{'seq'},"\t",$nomod_pep,"\t",$mod,"\t";
#				my $peptide_mass = 0;				
				
#				$peptide_mass = $massutil->get_peptide_mass($nomod_pep,$mod,"N");
				
	#			print $prec_mass,"\t",$peptide_mass,"\n";
	
########### In version2.0.1, we merged all individual hashes into one MatchedMass hash 	
				$MatchedMass{$PeptideID}{$mod}{'proteinID'} = $peptidehash->{$PeptideID}->{'proteinID'};
				$MatchedMass{$PeptideID}{$mod}{'origseq'} = $peptidehash->{$PeptideID}->{'seq'};
				$MatchedMass{$PeptideID}{$mod}{'seq'} = $seq_array[1];
				$MatchedMass{$PeptideID}{$mod}{'nomod_seq'} = $nomod_pep;
				$MatchedMass{$PeptideID}{$mod}{'theoretical_MH'} = $prec_mass;	

###################################################				
#				$proteinid{$PeptideID}{$mod} = $peptidehash->{$PeptideID}->{'proteinID'};			
#				$tagGroups_orig{$PeptideID}{$mod}=$peptidehash->{$PeptideID}->{'seq'};
#				$tagGroups{$PeptideID}{$mod}=$nomod_pep;
########### in version 2: we used $prec_mass as the theorectical mass even it is a rough value ############				
#				$threotical_MH{$PeptideID}{$mod}=$peptide_mass;
#				$threotical_MH{$PeptideID}{$mod} = $prec_mass;
			}			
		}
	}
	$self->{'matched_prec_peptide_num'} = scalar keys %MatchedMass;
	return \%MatchedMass;
}
	
############################ modification part has been removed because it is added to database ##################	
# if there are any dynamic modifications 
#	my $parameter = new Spiders::Params();
#	my $math = new Spiders::MathUtils();
	
#	my ($dynamic_modif, $largest_modif) = $parameter->get_dynamic_modifications($params);
#	my %dynamic_modif_comb = $math->combinations(\%$dynamic_modif,$max_modif);
#	foreach $modif_comb (keys %dynamic_modif_comb)
#	{
#		my $prec_start = int(($tags->{'precMass'}+$dynamic_modif_comb{$modif_comb}-$error)*1000);
#		my $prec_end = int(($tags->{'precMass'}+$dynamic_modif_comb{$modif_comb}+$error)*1000);
#		$self->search_core(\%proteinid,\%tagGroups,\%threotical_MH,$prec_start,$prec_end);		
#	}
#	my $mod_num = scalar keys %$dynamic_modif;	
######## get how many possible combinations ################################################################					
					
=head					
		my $prec_start1 = int(($tags->{'precMass'}-$max_modif * $largest_modif)*1000);
		my $prec_end1 = int(($tags->{'precMass'})*1000);	
		for(my $prec_mass = $prec_start1; $prec_mass <= $prec_end1; $prec_mass++)
		{
			my $massId = $index->readMassIndex($mass_index_file,$prec_mass);
			if($massId!=-1){
				my $peptidehash = $index->getAssociatedPeptides($massId,$peptide_index_file);
				foreach $PeptideID (keys %{$peptidehash})
				{		
					next if (!($self->matchtwoseq($peptidehash->{$PeptideID}->{'seq'},$tags->{'tagSeq'})));
					my @modif_com;
					my $length=length($peptidehash->{$PeptideID}->{'seq'});
					my $mod = ":" x $length;

					push (@modif_com,$mod);				

	###### count how many modif amino acid in the peptide		
					foreach my $modif_amino (keys %$modif)
					{
########### How many amino acid matched??????????????

						while($peptidehash->{$PeptideID}->{'seq'} =~ /$modif_amino/g)
						{		
							my @modif_com_old = @modif_com;
							for(my $i=0;$i<=$#modif_com_old;$i++)
							{

								substr($modif_com_old[$i],$-[0],1)= ":" . $modif_amino ;						
								push(@modif_com,$modif_com_old[$i]);
							}						
						}
					}
				}
		####### remove the unmodif	
				shift @modif_com;
	#		next if (scalar ($#modif_com) ==0);
	######### if there is any modif amino acids, then what the combination? #####
				for(my $j=0;$j<=$#modif_com;$j++)
				{
					my $peptide_mass = $massutil->get_peptide_mass($peptidehash->{$PeptideID}->{'seq'},$modif_com[$j],"N");
					if (abs($peptide_mass - $tags->{'precMass'}) <= $error*1000){				
						$proteinid{$PeptideID}{$mod} = $peptidehash->{$PeptideID}->{"proteinID"};	
						$tagGroups{$PeptideID}{$modif_com[$j]}=$peptidehash->{$PeptideID}->{'seq'};
						$threotical_MH{$PeptideID}{$modif_com[$j]}=$peptide_mass;
					}				
				}
			}
		
		}
	}
=cut	
########## version #################################################################


=head	
	foreach my $mass (keys %{$masshash})
	{
		next if (!defined($masshash->{$mass}));
		
		
	#	my $mw = substr($database->[$i]->{'id'},0,index($database->[$i]->{'id'},'|'));	
## when matching precursor mass and theoritcal mass, the formula is Mass(prec) = Mass(theo) + H
		if (abs($mass - $tags->{'precMass'}*1000) <= $error*1000){	
## record the mass error

			foreach my $PeptideID (@{$masshash->{$mass}})
			{
				
				my $length=length($peptidehash->{$PeptideID}->{'seq'});
				my $mod = ":" x $length;		
				$tagGroups{$PeptideID}{$mod}=$peptidehash->{$PeptideID}->{'seq'};
				$precursor{$PeptideID}{$mod}=$mass;
			}
		}
=cut		
######## MS1 isotopic peaks
=head MS1 deisotoping performed in decharge.pm		
		elsif(abs($mass + $H*1000 - $tags->{'precMass'}*1000) <= $error*1000)
		{
## record the mass error
			foreach my $PeptideID (@{$masshash->{$mass}})
			{
				my $length=length($peptidehash->{$PeptideID}->{'seq'});
				my $mod = ":" x $length;
				$tagGroups{$PeptideID}{$mod}=$peptidehash->{$PeptideID}->{'seq'};
				$precursor{$PeptideID}{$mod}=$mass;
			}
		}
=cut
=head		
		elsif($mass < $tags->{'precMass'}*1000 and ($tags->{'precMass'}*1000 < ($mass +($max_modif * $largest_modif*1000))))  
		{
### it will dramatically reduce  the run time by searching tags

			foreach my $PeptideID (@{$masshash->{$mass}})
			{	
				next if (!($self->matchtwoseq($peptidehash->{$PeptideID}->{'seq'},$tags->{'tagSeq'})));
				my @modif_com;
				my $length=length($peptidehash->{$PeptideID}->{'seq'});
				my $mod = ":" x $length;

				push (@modif_com,$mod);				

	###### count how many modif amino acid in the peptide		
				foreach my $modif_amino (keys %$modif)
				{
########### How many amino acid matched??????????????

					while($peptidehash->{$PeptideID}->{'seq'} =~ /$modif_amino/g)
					{		
						my @modif_com_old = @modif_com;
						for(my $i=0;$i<=$#modif_com_old;$i++)
						{

							substr($modif_com_old[$i],$-[0],1)= ":" . $modif_amino ;						
							push(@modif_com,$modif_com_old[$i]);
						}						
					}
				}
			}
		####### remove the unmodif	
			shift @modif_com;
	#		next if (scalar ($#modif_com) ==0);
	######### if there is any modif amino acids, then what the combination? #####
			for(my $j=0;$j<=$#modif_com;$j++)
			{
				my $peptide_mass = $massutil->get_peptide_mass($peptidehash->{$PeptideID}->{'seq'},$modif_com[$j],"N");
				if (abs($peptide_mass - $tags->{'precMass'}) <= $error){				
					$tagGroups{$PeptideID}{$modif_com[$j]}=$peptidehash->{$PeptideID}->{'seq'};
					$precursor{$PeptideID}{$modif_com[$j]}=$peptide_mass;
				}				
			}
		}
	}
=cut

####### In version 2.0.1, we change "findTags" function into two independent functions: SearchMass and SearchTag. 
####### The advantage of sperating into two functions is only search once for mass no matter how many tags will be used for searching, where is especially useful for low resolution.
 

sub SearchTag
{	
	my ($self,$MatchedMassRef) = @_;

	my $matchResults;
	
	my $index = new Spiders::BuildIndex();
	my $database = $self->get_database();	
	my $protein_index_file = $database . ".prdx";	
	
	my $massutil = new Spiders::MassUtils();
	my $params = $self->get_parameter();	
	$massutil->set_parameter($params);
	my $tags = $self->get_tag();
	my $frag_tolerance = $self->get_frag_tolerance();
	
	
	my %MatchedMass = %$MatchedMassRef;
	
	foreach my $pep (keys %MatchedMass)
	{
		foreach $mod (keys %{$MatchedMass{$pep}})
		{
			my $pepseq = $MatchedMass{$pep}{$mod}{'nomod_seq'};
			my $pepseq_core = $MatchedMass{$pep}{$mod}{'nomod_seq'};
#			my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq,$mod,$tags,$error);
			my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq_core,$mod,$tags,$frag_tolerance);
#			print $pepseq_core,"\t",$sideMassMatchResult->{'lSide_match'},"\t",$sideMassMatchResult->{'rSide_match'},"\t",$tags->{'tagSeq'},"\n";
			if ($sideMassMatchResult->{'lSide_match'} || $sideMassMatchResult->{'rSide_match'})
			{	
				my $protein_id = $MatchedMass{$pep}{$mod}{'proteinID'};
				
				my ($proteinName,$proteinDesc) = $index->getProtein($protein_id,$protein_index_file);
							
				$matchResults->{$pep}->{$mod}->{'scanNum'} = $tags->{'scanNum'};
				$matchResults->{$pep}->{$mod}->{'precMass'} = $tags->{'precMass'};
				$matchResults->{$pep}->{$mod}->{'theoretical_MH'} = $MatchedMass{$pep}{$mod}{'theoretical_MH'};
				$matchResults->{$pep}->{$mod}->{'tagSeq'} = $tags->{'tagSeq'};			
				$matchResults->{$pep}->{$mod}->{'sideMass'} = $tags->{'sideMass'};			
				$matchResults->{$pep}->{$mod}->{'pepseq'} = $pepseq;
				$matchResults->{$pep}->{$mod}->{'proteinid'} = $proteinName;

												
				$matchResults->{$pep}->{$mod}->{'tagSeqMatch'} = $tagSeqMatchResult;
				$matchResults->{$pep}->{$mod}->{'lsideMassMatch'} = $sideMassMatchResult->{'lSide_match'};
				$matchResults->{$pep}->{$mod}->{'rsideMassMatch'} = $sideMassMatchResult->{'rSide_match'};
				$matchResults->{$pep}->{$mod}->{'lthreot_Mass'} = sprintf("%.4f",$sideMassMatchResult->{'lthreot_Mass'});
				$matchResults->{$pep}->{$mod}->{'rthreot_Mass'} = sprintf("%.4f",$sideMassMatchResult->{'rthreot_Mass'});

				$matchResults->{$pep}->{$mod}->{'matchPeptide'} = $pepseq;

				$matchResults->{$pep}->{$mod}->{'matchPeptide_orig'} = $MatchedMass{$pep}{$mod}{'origseq'};

				$matchResults->{$pep}->{$mod}->{'rank_p'} = $tags->{'rankP'};
				$matchResults->{$pep}->{$mod}->{'hyper_p'} = $tags->{'hyperP'};
			}
			elsif($tagSeqMatchResult eq '0')
			{
				my $combination_mass = $massutil->get_AA_combination();
				foreach my $aa (keys %$combination_mass)
				{
					if($tags->{'tagSeq'} =~ /$aa/)
					{
						my $replace_aa = $combination_mass->{$aa};
						foreach $replace_aa (@$replace_aa)
						{
							my %tags_temp = %$tags;
							$tags_temp{'tagSeq'} =~ s/$aa/$replace_aa/;
							my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq_core,$mod,\%tags_temp,$frag_tolerance);
							if ($sideMassMatchResult->{'lSide_match'} || $sideMassMatchResult->{'rSide_match'})
							{			
							
								my ($proteinName,$proteinDesc) = $index->getProtein($MatchedMass{$pep}{$mod}{'proteinID'},$protein_index_file);							
								
								$matchResults->{$pep}->{$mod}->{'scanNum'} = $tags->{'scanNum'};
								$matchResults->{$pep}->{$mod}->{'precMass'} = $tags->{'precMass'};
								$matchResults->{$pep}->{$mod}->{'theoretical_MH'} = $MatchedMass{$pep}{$mod}{'theoretical_MH'};
						####### need to fix blank bug if use the updated tag seq ################		
								# $matchResults->{$pep}->{$mod}->{'tagSeq'} = $tags_temp->{'tagSeq'};
								 $matchResults->{$pep}->{$mod}->{'tagSeq'} = $tags->{'tagSeq'};
								
								$matchResults->{$pep}->{$mod}->{'sideMass'} = $tags->{'sideMass'};			
								$matchResults->{$pep}->{$mod}->{'pepseq'} = $pepseq;
								$matchResults->{$pep}->{$mod}->{'proteinid'} = $proteinName;	
								$matchResults->{$pep}->{$mod}->{'tagSeqMatch'} = $tagSeqMatchResult;
								$matchResults->{$pep}->{$mod}->{'lsideMassMatch'} = $sideMassMatchResult->{'lSide_match'};
								$matchResults->{$pep}->{$mod}->{'rsideMassMatch'} = $sideMassMatchResult->{'rSide_match'};
								$matchResults->{$pep}->{$mod}->{'lthreot_Mass'} = sprintf("%.4f",$sideMassMatchResult->{'lthreot_Mass'});
								$matchResults->{$pep}->{$mod}->{'rthreot_Mass'} = sprintf("%.4f",$sideMassMatchResult->{'rthreot_Mass'});

								$matchResults->{$pep}->{$mod}->{'matchPeptide'} = $pepseq;
								$matchResults->{$pep}->{$mod}->{'matchPeptide_orig'} = $MatchedMass{$pep}{$mod}{'origseq'};							
								$matchResults->{$pep}->{$mod}->{'rank_p'} = $tags->{'rankP'};
								$matchResults->{$pep}->{$mod}->{'hyper_p'} = $tags->{'hyperP'};
							}
						}						
					}
				}				
			}
		}
	}
	my $num_peptides = scalar keys %$matchResults;

	
	if($num_peptides > 0)
	{
		return $matchResults;
	}
	else
	{
		return 0;
	}
}





sub matchtwoseq
{
	my ($self,$seq1,$seq2)=@_;
	my $matched=0;
	my $pseq = $seq1;
	$pseq =~ s/I/J/g;
	$pseq =~ s/L/J/g;

	my $tseq = $seq2;
	$tseq =~ s/I/J/g;
	$tseq =~ s/L/J/g;
	my $rev_tseq=reverse($tseq);
	if($pseq=~/$tseq/)
	{
		$matched=1;
	}
	elsif($pseq=~/$rev_tseq/)
	{
		$matched=1;
	}
	else
	{
		$matched=0;
	}
	return $matched;
}

sub matchTagSeq{

	my ($self,  $pepseq, $mod, $tags, $error) = @_;


	my $pseq = $pepseq;
	$pseq =~ s/I/J/g;
	$pseq =~ s/L/J/g;

	my $tseq = $tags->{'tagSeq'};
	$tseq =~ s/I/J/g;
	$tseq =~ s/L/J/g;
	
	my $tagSeqMatch = 0;
	my $sideMassMatch ;
	my $pepsideMass = 0;
	my $whole = -1;
	my $lPart = -1;
	my $rPart = -1;
	my $removed_seq = "";

	while (1) {	
		$whole = index($pseq,$tseq,$whole+1);
		
		if($whole < 0)
		{
			$sideMassMatch->{'lSide_match'} = 0;
			$sideMassMatch->{'rSide_match'} = 0;
			last;			
		}
		$tagSeqMatch = 1;
		$sideMassMatch = $self->checkSideMass($whole,$tags,$pseq,$mod,$error,$removed_seq);
		last if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0);
	}
=head	###### ranking tag by p_value rather than using this removing left and right aa method
	if($sideMassMatch->{'lSide_match'}==0 and $sideMassMatch->{'rSide_match'}==0)
	{
		while (1) {	
			$lPart = (length($tseq) > 1) ? index($pseq, substr($tseq,0,length($tseq)-1),$lPart+1) : -1;
			last if($lPart < 0);
			$tagSeqMatch = 2;
#			$removed_seq = substr($tseq,length($tseq)-1,1);
			
			$sideMassMatch = $self->checkSideMass($lPart,$tags,$pseq,$mod,$error,$removed_seq);
			last if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0);
		}	
	}
	if($sideMassMatch->{'lSide_match'}==0 and $sideMassMatch->{'rSide_match'}==0)
	{
		while (1) {	
			$rPart = (length($tseq) > 1) ? index($pseq, substr($tseq,1,length($tseq)-1),$rPart+1) : -1;
			last if($rPart < 0);
			$tagSeqMatch = 3;
			$removed_seq = substr($tseq,0,1);
			$sideMassMatch = $self->checkSideMass($rPart,$tags,$pseq,$mod,$error,$removed_seq);
			last if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0);
		}	
	}
=cut
	$whole = -1;
	if($sideMassMatch->{'lSide_match'}==0 and $sideMassMatch->{'rSide_match'}==0)
	{
		my $rev_tseq = reverse ($tseq);

		while (1) {	
			$whole = index($pseq,$rev_tseq,$whole+1);
		

			last if($whole < 0);
			$tagSeqMatch = 1;
	
			$sideMassMatch = $self->checkSideMass_rev($whole,$tags,$pseq,$mod,$error,$removed_seq);
			last if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0);	
		}
	}
=head	###### ranking tag by p_value rather than using this removing left and right aa method
	if($sideMassMatch->{'lSide_match'}==0 and $sideMassMatch->{'rSide_match'}==0)
	{
		my $rev_tseq = reverse ($tseq);
		while (1) {	
			$lPart = (length($rev_tseq) > 1) ? index($pseq, substr($rev_tseq,0,length($rev_tseq)-1),$lPart+1) : -1;
			last if($lPart < 0);
			$tagSeqMatch = 2;
#				$removed_seq = substr($tseq,length($tseq)-1,1);
				
			$sideMassMatch = $self->checkSideMass_rev($lPart,$tags,$pseq,$mod,$error,$removed_seq);
			last if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0);
		}	
	}
	if($sideMassMatch->{'lSide_match'}==0 and $sideMassMatch->{'rSide_match'}==0)
	{
		my $rev_tseq = reverse ($tseq);
		while (1) {	
			$rPart = (length($rev_tseq) > 1) ? index($pseq, substr($rev_tseq,1,length($rev_tseq)-1),$rPart+1) : -1;		
				last if($rPart < 0);
				$tagSeqMatch = 3;
				$removed_seq = substr($tseq,0,1);
				$sideMassMatch = $self->checkSideMass_rev($rPart,$tags,$pseq,$mod,$error,$removed_seq);
				last if($sideMassMatch->{'lSide_match'} != 0 || $sideMassMatch->{'rSide_match'} != 0);	
		}	
		
	}
=cut	
	if(!defined($sideMassMatch->{'lSide_match'}))
	{
		$sideMassMatch->{'lSide_match'}=0;
	}
	if(!defined($sideMassMatch->{'rSide_match'}))
	{
		$sideMassMatch->{'rSide_match'}=0;
	}
	
	return ($tagSeqMatch,$sideMassMatch);
	
}

sub checkSideMass{

	my ($self, $ind, $tag, $pseq, $mod, $error,$removed_seq) = @_;
	my $massutil = new Spiders::MassUtils();
	
	my $parameter = $self->get_parameter();
	$massutil->set_parameter($parameter);
	
######## N-term mass has to add H2O + H
	my $result;
	my $pepSideMass;

	my $N_term_Mass = 19.017806;
	my $lSide = substr($pseq, 0, $ind);
	my $rSide = substr($pseq, $ind);
	my $lSide_mod = substr($mod,0,$ind);
	my $rSide_mod = substr($mod, $ind);
	if($removed_seq ne "")
	{
	##### need check ??????##
	#	$tag->{'sideMass'} = $tag->{'sideMass'} - $massutil->get_peptide_mass($removed_seq,"","C") + 1.007825;
	}
####### because the tag was reverse, see tag.pm module, so the N and C terminus was reversed

	$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"C");

	$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"N");

#	print $pseq,"\t",$lSide,"\t",$rSide,"\t",$tag->{'sideMass'},"\t",$result->{'lthreot_Mass'},"\n"; 
	
### commented by xusheng; only error without multiple with $tag->{'sideMass'}
	
	if (abs($tag->{'sideMass'} - $result->{'lthreot_Mass'} ) <= $error )
	{
		$result->{'lSide_match'} =  1;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};
	}
	else
	{
#		$result->{'lthreot_Mass'} = $self->calMW($lSide);

		$result->{'lSide_match'} = 0;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};	
	}

	if (abs($result->{'rthreot_Mass'} - $tag->{'sideMass'}) <= $error ){
		$result->{'rSide_match'} = 1 ;

		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}
	else{
		$result->{'rSide_match'} = 0 ;
		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}
=head	
######## This is not the right case ##########
	if($result->{'lSide_match'} eq '0' and $result->{'rSide_match'} eq '0')
	{
		$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"N");
		$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"C");	
		if (abs($tag->{'sideMass'} - $result->{'lthreot_Mass'} ) <= $error ){

			$result->{'lSide_match'} =  1;
			$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};
		}
		else
		{
	#		$result->{'lthreot_Mass'} = $self->calMW($lSide);
			$result->{'lSide_match'} = 0;
			$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};	
		}

		if (abs($result->{'rthreot_Mass'} - $tag->{'sideMass'}) <= $error ){
			$result->{'rSide_match'} = 1 ;

			$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
		}
		else{
			$result->{'rSide_match'} = 0 ;
			$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
		}		
	}
=cut 	
	
	return $result;
}

sub checkSideMass_rev{

	my ($self, $ind, $tag, $pseq, $mod, $error,$removed_seq) = @_;
	my $massutil = new Spiders::MassUtils();
	
	my $parameter = $self->get_parameter();
	$massutil->set_parameter($parameter);
	
######## N-term mass has to add H2O + H
	my $result;
	my $pepSideMass;

## to match the side Mass, needs to add the length of tag because N and C terminus shift;  	
	$ind += length($tag->{'tagSeq'});
#####	
	my $lSide = substr($pseq, 0, $ind);
	my $rSide = substr($pseq, $ind);

	
	my $lSide_mod = substr($mod,0,$ind);
	my $rSide_mod = substr($mod, $ind);
	if($removed_seq ne "")
	{
#		$tag->{'sideMass'} = $tag->{'sideMass'} - $massutil->get_peptide_mass($removed_seq,"","C") + 1.007825;
	}
####### because the tag was reverse, see tag.pm module, so the N and C terminus was reversed
########### ????????????????????????????????
	$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"C") - 19.017806 + 1.007825;

	$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"N") + 19.017806 - 1.007825;
########## ??????????????????????????????????????

	
### commented by xusheng; only error without multiple with $tag->{'sideMass'}
#	print $pseq,"\t",$lSide,"\t",$rSide,"\t",$tag->{'sideMass'},"\t",$result->{'lthreot_Mass'},"\t",$result->{'rthreot_Mass'},"uuuuuuuu\n";	
	if (abs($tag->{'sideMass'} - $result->{'lthreot_Mass'} ) <= $error ){

		$result->{'lSide_match'} =  1;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};
	}
	else
	{
#		$result->{'lthreot_Mass'} = $self->calMW($lSide);
		$result->{'lSide_match'} = 0;
		$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};	
	}

	if (abs($result->{'rthreot_Mass'} - $tag->{'sideMass'}) <= $error ){
		$result->{'rSide_match'} = 1 ;

		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}
	else{
		$result->{'rSide_match'} = 0 ;
		$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
	}
=head	
######## This is not the right case ##########	
	if($result->{'lSide_match'} eq '0' and $result->{'rSide_match'} eq '0')
	{
		$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"N");
		$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"C");	
		if (abs($tag->{'sideMass'} - $result->{'lthreot_Mass'} ) <= $error ){

			$result->{'lSide_match'} =  1;
			$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};
		}
		else
		{
	#		$result->{'lthreot_Mass'} = $self->calMW($lSide);
			$result->{'lSide_match'} = 0;
			$result->{'lSide_masserror'} = $result->{'lthreot_Mass'} - $tag->{'sideMass'};	
		}

		if (abs($result->{'rthreot_Mass'} - $tag->{'sideMass'}) <= $error ){
			$result->{'rSide_match'} = 1 ;

			$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
		}
		else{
			$result->{'rSide_match'} = 0 ;
			$result->{'rSide_masserror'} = $result->{'rthreot_Mass'} - $tag->{'sideMass'};
		}	
	}
=cut
	return $result;
}

=head
sub calMW{

	my ($self, $sq) = @_;
	my $parameter = $self->get_parameter();
	
	if (length($sq) ==0){
		return 0;
	}else{
		my $aaMW = {'J' => 113.08406,
					'I' => 113.08406,
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

######### Get modification 

		
		my $h = {};
		my @aas = split(//,$sq);
		foreach (@aas){
			$h->{$_}++;
		}
		my $mw = 0;
		foreach $aa (keys %$h){
			if(!defined($aaMW->{$aa}))
			{
	#			print "\ryour sequences contain unkown amino acid: $aa!";
				next;
			}
			$mw += $h->{$aa} * $aaMW->{$aa};
		}
		return $mw;
	}
}
=cut

############## Version 10.09 ######################
## Add the function allowing to search without Tag matches 

sub SearchWithoutTag
{	
	my ($self,$MatchedMassRef) = @_;

	my $matchResults;
	
	my $index = new Spiders::BuildIndex();
	my $database = $self->get_database();	
	my $protein_index_file = $database . ".prdx";	
	
	my $massutil = new Spiders::MassUtils();
	my $params = $self->get_parameter();	
	$massutil->set_parameter($params);
	my $tags = $self->get_tag();
	my $frag_tolerance = $self->get_frag_tolerance();
	
	
	my %MatchedMass = %$MatchedMassRef;
	
	foreach my $pep (keys %MatchedMass)
	{
		foreach $mod (keys %{$MatchedMass{$pep}})
		{
			my $pepseq = $MatchedMass{$pep}{$mod}{'nomod_seq'};


			my $protein_id = $MatchedMass{$pep}{$mod}{'proteinID'};
				
			my ($proteinName,$proteinDesc) = $index->getProtein($protein_id,$protein_index_file);
							
			$matchResults->{$pep}->{$mod}->{'scanNum'} = $tags->{'scanNum'};
			$matchResults->{$pep}->{$mod}->{'precMass'} = $tags->{'precMass'};
			$matchResults->{$pep}->{$mod}->{'theoretical_MH'} = $MatchedMass{$pep}{$mod}{'theoretical_MH'};
			$matchResults->{$pep}->{$mod}->{'tagSeq'} = "N/A";			
			$matchResults->{$pep}->{$mod}->{'sideMass'} = "";			
			$matchResults->{$pep}->{$mod}->{'pepseq'} = $pepseq;
			$matchResults->{$pep}->{$mod}->{'proteinid'} = $proteinName;

												
			$matchResults->{$pep}->{$mod}->{'tagSeqMatch'} = "";
			$matchResults->{$pep}->{$mod}->{'lsideMassMatch'} = "";
			$matchResults->{$pep}->{$mod}->{'rsideMassMatch'} = "";
			$matchResults->{$pep}->{$mod}->{'lthreot_Mass'} = "";
			$matchResults->{$pep}->{$mod}->{'rthreot_Mass'} = "";

			$matchResults->{$pep}->{$mod}->{'matchPeptide'} = $pepseq;

			$matchResults->{$pep}->{$mod}->{'matchPeptide_orig'} = $MatchedMass{$pep}{$mod}{'origseq'};
			$matchResults->{$pep}->{$mod}->{'rank_p'} = "";
			$matchResults->{$pep}->{$mod}->{'hyper_p'} = "";
		}
	}
	my $num_peptides = scalar keys %$matchResults;
	if($num_peptides > 0)
	{
		return $matchResults;
	}
	else
	{
		return 0;
	}	
}

sub generate_peptide_theoritical_mass
{
	my ($self,$peptide,$modif)=@_;
	my $parameter = $self->get_parameter();
	my $massutils = new Spiders::MassUtils();
	$massutils->set_parameter($parameter);
######## If there is static modification ################
# no need to change the $modif variable because it has been set in AA_mass
=head
	my $param = new Spiders::Params();
	my $static_modif = $param->get_static_modification($parameter);

	
	foreach my $static_mod (keys %$static_modif)
	{
		my $pos = -1;
		while (1) {
			$pos = index($peptide, $static_mod, $pos + 1);
			last if $pos < 0;

			$modif = $modif . "-";
			my @modif_array = split(/:/,$modif);
			$modif_array[$pos+1] .= $static_mod;
			$modif = join(":",@modif_array);
			$modif =~ s/\-//g;	
		}
	}
=cut	
##########################################################	
##########Version 1.08########################	
# add the ion series function 	
	my $ion_series = $parameter->{'ion_series'};
	my @ions_array = split(/\s+/,$ion_series);
	if(scalar(@ions_array) != 9)
	{
		print "please put the right ion_series parameter!\n";
		exit;
	}
	my @ion_series_name = ('a', 'b', 'c', 'd', 'v', 'w', 'x', 'y', 'z','a++', 'b++', 'c++', 'd++', 'v++', 'w++', 'x++', 'y++', 'z++');
	my @ion_seies_used;
	for (my $i=0;$i<$#ions_array;$i++)
	{
		if($ions_array[$i]==1)
		{
			push (@ion_seies_used,$ion_series_name[$i]);
		}
	}
# add the charge state 
	my $charge = $self->get_precCharge();
	if($charge>=2)
	{
		for (my $i=0;$i<$#ions_array;$i++)
		{	
			if($ions_array[$i]==1)
			{
				push (@ion_seies_used,$ion_series_name[$i+9]);
			}
		}
	}
	
	$massutils->getFragmentMasses(pept=>$peptide, modif=>$modif, fragTypes=>[@ion_seies_used], spectrum=>\%spectrum);
#	$massutils->getFragmentMasses(pept=>$peptide, modif=>$modif, fragTypes=>['b','y','b++','y++','immo'], spectrum=>\%spectrum);
####################end ####################################

	my (@aa, @bb);
	
	foreach (%{$spectrum{'mass'}{'term'}}){

		push (@aa, @{$spectrum{'mass'}{'term'}{$_}});
		
	}
	foreach (%{$spectrum{'mass'}{'intern'}}){
		@bb =  @{$_->{'masses'}};
	}

	my @big = (@aa,@bb);
=head
	foreach (@aa)
	{
		print $_,"aaaaaaa\n";
	}
	foreach (@bb)
	{
		print $_,"bbbbbbb\n";
	}
=cut	
	return (\@big);
	
}

sub compare_theoritical_experiment
{
	my ($self,$exp_mz,$theor_mz) = @_;
	my $matched = 0;
	my $tolerance = $self->get_frag_tolerance();

	foreach my $exp_mz_value (@$exp_mz)
	{
		foreach my $theor_mz_value (@$theor_mz)
		{

			if(($exp_mz_value+$tolerance)>$theor_mz_value and ($exp_mz_value-$tolerance)<$theor_mz_value)
			{
#				print $theor_mz_value,"\t",$exp_mz_value,"ffffffffff\n";
				$matched++;
			}
#			print $theor_mz_value,"xxxxxxxxxxxx\n";		
		}
	}
	return $matched;
}

sub set_exp_mz
{
	my ($self,$exp_mz)=@_;
	$self->{'_exp_mz'} = $exp_mz;
}

sub get_exp_mz
{
	my $self=shift;
	return $self->{'_exp_mz'};
}

=head
sub get_peptide_pvalue
{
	my ($self,$peptide,$modif)=@_;
	my $theor_mz = $self->generate_peptide_theoritical_mass($peptide,$modif);
	my $exp_mz = $self->get_exp_mz();
	my $matched = $self->compare_theoritical_experiment($exp_mz,$theor_mz);
	
	my $mass_tolerance = $self->get_mass_tolerance();
	
	my $total_location = ($#$theor_mz+1)/$mass_tolerance;
	
	my $theor_mass_num = scalar (@$theor_mz);
	my $exp_mass_num = scalar (@$exp_mz);
	
	my $hyper = new Spiders::Hypergeometric();
	my $peptide_pvalue=$hyper->Hypergeometric($total_location,$theor_mass_num,$exp_mass_num,$matched);

	my $log_peptide_pvalue = $peptide_pvalue>0 ? sprintf("%.4f",-log($peptide_pvalue)) : 0;

############ Version 1.08 ########################### 	
# add return theoretical ions and matched ions
	$self->{'Theor_ion_num'} = $#$theor_mz + 1; 
	$self->{'Matched_ion_num'} = $matched; 
################### end ##############################	

	return $log_peptide_pvalue;
}
=cut

#### Version 10.11 in order to increase the speed #####################

sub get_peptide_matched
{
	my ($self,$peptide,$modif)=@_;
	my $theor_mz = $self->generate_peptide_theoritical_mass($peptide,$modif);
	my $exp_mz = $self->get_exp_mz();
#	print $peptide,"\n";
	my $matched = $self->compare_theoritical_experiment($exp_mz,$theor_mz);

############ Version 1.08 ########################### 	
# add return theoretical ions and matched ions
	$self->{'Theor_ion_num'} = $#$theor_mz + 1; 
	$self->{'Matched_ion_num'} = $matched; 
################### end ##############################	

	return $matched;
}

sub get_peptide_pvalue
{
	my ($self,$matched,$total)=@_;
###### Version 1.12 fix a bug: using frag tolerance instead of mass tolerance	
	my $mass_tolerance = $self->get_frag_tolerance();
	
	my $total_location = $total/$mass_tolerance;
	my $exp_mz = $self->get_exp_mz();
	
	
#	foreach (@$exp_mz)
#	{
#		print $_,"\n";
#	}


	my $exp_mass_num = scalar (@$exp_mz);
	
	my $hyper = new Spiders::Hypergeometric();
#	print $total_location,"\t",$total,"\t",$exp_mass_num,"\t",$matched,"?????????????????\n";
	my $peptide_pvalue=$hyper->Hypergeometric($total_location,$total,$exp_mass_num,$matched);

	my $log_peptide_pvalue = $peptide_pvalue>0 ? sprintf("%.4f",-log($peptide_pvalue)/log(10)) : 0;

################### end ##############################	

	return $log_peptide_pvalue;
}



sub mergeresults
{
	my ($self,$result1,$result2,$tag) = @_;
	foreach (keys %$result2)
	{
		$result1->{$tag}->{$_} = $result2->{$_};
	}
	return $result1;
}

sub writingHeader
{
	my ($self,$outfile) = @_;
	my $parameter = $self->get_parameter();
## multiple tag #######	
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";

	print OUTPUT "\nspTag version 2.01 (c) 2012\n\n";
	print OUTPUT "St Jude Children's Research Hospital, X Wang/J Peng\n";
	my $dtafile = $self->get_scanNum();
	$dtafile =~ s/tag/dta/;
	print OUTPUT "DTA file = ", $dtafile,"\n";
	print OUTPUT "Database =", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass =",$self->get_precMass(),"\n\n";
	close(OUTPUT);
	
}

sub WriteResults
{
	my ($self,$outfile,$matches,$tag_rank_num,$prec_peak_int) = @_;
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";
	my $parameter = $self->get_parameter();
	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");	
## multiple tag #######	
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";


	print OUTPUT "\nJUMP version 2.01 (c) 2012\n\n";
	print OUTPUT "St Jude Children's Research Hospital, X Wang/J Peng\n";
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	printf OUTPUT "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon,$mday,$hour,$min,$sec;	
	my $dtafile = $self->get_scanNum();
	$dtafile =~ s/tag/dta/;
	print OUTPUT "DTA file = ", $dtafile,"\n";
	
	print OUTPUT "Database = ", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass = ",$self->get_precMass(),"   ", "Percentage of precursor peak intensity = ", sprintf("%0.2f",$prec_peak_int*100),"\%\n";
	print OUTPUT "Precursor matched peptides = ", $self->{'matched_prec_peptide_num'},"\, ","Peptides mass tolerance = ", $self->get_mass_tolerance(), "  Peptides mass tolerance units = $parameter->{'peptide_tolerance_units'}\n";
	print OUTPUT "Tag number = ",$self->get_tag_number(),"\n";
	
	print OUTPUT "ion series ABCDVWXYZ: ", $parameter->{'ion_series'},"\n";

###### output modification information ###############################
	my $i=0;
	foreach (keys %$parameter)
	{
		if($_ =~/add_/ and $parameter->{$_} != 0)
		{
			print OUTPUT $_,"=",$parameter->{$_},"  ";
		}

		if($_ =~/dynamic_/ and $parameter->{$_} != 0)
		{
			print OUTPUT $_,$mod_symbol[$i],"=",$parameter->{$_},"  ";
			$i++;
		}
	}
	
	print OUTPUT "\n\n";
#	print OUTPUT "Tag = ", $self->get_tagSeq(),"\n";
#	my $tag_e_value = sprintf("%.2f",$self->get_tagPvalue());
#	print OUTPUT "Tag E_value = ", $tag_e_value,"\n","Tag Rank Number = ",$tag_rank_num,"\n","Tag Side Mass = ",$self->get_sideMass(),"\n"; 
	my $i=0;
	my %sort_results;
	my %comb_results;
	my %peptide_matched_ratio;
	my %ions;
	
	foreach my $tag (keys %{$matches})
	{
		foreach my $pep (keys %{$matches->{$tag}})
		{	
			my $j=1;
			foreach my $mod (keys %{$matches->{$tag}->{$pep}})
			{	
#					print $tag,"\t",$pep,"\t",$mod,"\n"; 
	#		print $matches->{$pep}->{'tagSeqMatch'},"\n";
					next if ($matches->{$tag}->{$pep}->{$mod}->{'tagSeqMatch'} < 1);
		#		print $matches->{$pep}->{$mod}->{'lsideMassMatch'},"\t",$matches->{$pep}->{$mod}->{'rsideMassMatch'},"\n";
	#			if ((abs($matches->{$pep}->{$mod}->{'lthreot_Mass'} - $self->get_sideMass()) < 0.02) || (abs($matches->{$pep}->{$mod}->{'rthreot_Mass'} - $self->get_sideMass()) < 0.02))
	#			{
					my $peptideseq = $matches->{$tag}->{$pep}->{$mod}->{'pepseq'};
#					if(!$peptide_pvalue{$pep})
					if(!$peptide_matched_ratio{$pep})
					{					
#						my $pep_pvalue = $self->get_peptide_pvalue($peptideseq,$mod) + $matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}};
####### Version 10.10; increase search speed 
#						my $pep_pvalue = $self->get_peptide_pvalue($peptideseq,$mod);
						my $matched = $self->get_peptide_matched($peptideseq,$mod);
						next if ($matched==0);
						my $pep_matched_ratio = $self->{'Matched_ion_num'} / $self->{'Theor_ion_num'};

#						$peptide_pvalue{$pep} = $pep_pvalue;
						$peptide_matched_ratio{$pep} =  $pep_matched_ratio;	
						$ions{$pep} = "$self->{'Matched_ion_num'}\/$self->{'Theor_ion_num'}";
					}
	# version 1.02 $matches->{$tag}->{$pep}->{$mod}->{'rank_p'},
	# version 1.03 $matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}}
					my $results = sprintf("%-0.3f    %-0.4f    %-0.4f    %-0.2f     %-10s   %-0.2f    %-0.4f  %-10s  %-30s    %-20s", 
					$matches->{$tag}->{$pep}->{$mod}->{'theoretical_MH'}/1000,$matches->{$tag}->{$pep}->{$mod}->{'lthreot_Mass'},$matches->{$tag}->{$pep}->{$mod}->{'rthreot_Mass'},
					$peptide_matched_ratio{$pep},$matches->{$tag}->{$pep}->{$mod}->{'tagSeq'},$matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}},$matches->{$tag}->{$pep}->{$mod}->{'sideMass'},$ions{$pep},$matches->{$tag}->{$pep}->{$mod}->{'proteinid'},
					$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'});
					$sort_results{$peptide_matched_ratio{$pep}}{$tag}{$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'}}{$matches->{$tag}->{$pep}->{$mod}->{'proteinid'}}=$results;


#					$comb_results{$pep}{$tag}{'peptide_E_value'} =  $peptide_pvalue{$pep};
#					$comb_results{$pep}{$tag}{'peptide_seq'} =  $matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'};
#					$comb_results{$pep}{$tag}{'protein_id'} = $matches->{$tag}->{$pep}->{$mod}->{'proteinid'};
#					$comb_results{$pep}{$tag}{'tag_seq'} = $matches->{$tag}->{$pep}->{$mod}->{'tagSeq'};
#					$comb_results{$pep}{$tag}{'tag_E_value'} = $matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}};
#			}
			}
		}
	}


	
=head	
	my %selected;
	foreach my $peptide (keys %comb_results)
	{
		my $first_time = 1;	
		foreach my $tag (keys %{$comb_results{$peptide}})
		{
			if($first_time)
			{
				$selected{$peptide}{'pep_E_value'} = $comb_results{$pep}{$tag}{'pep_E_value'} + $comb_results{$pep}{$tag}{'tag_E_value'};
				$selected{$peptide}{'tag_seq'} = $comb_results{$pep}{$tag}{'tag_seq'};
				$selected{$peptide}{'peptide_seq'} = $comb_results{$pep}{$tag}{'peptide_seq'};
				$first_time = 0;
			}
			if($selected{$peptide}{'tag_seq'} =~ /$comb_results{$pep}{$tag}{'tag_seq'}/ || $comb_results{$peptide}{$tag}{'tag_seq'} =~ /$selected{$selected_pep}{'tag_seq'}/)
			{
				$selected{$peptide}{'pep_E_value'} +=  $comb_results{$peptide}{$tag}{'tag_E_value'} * 0.5;
			}
			else
			{
				$selected{$peptide}{'pep_E_value'} +=  $comb_results{$peptide}{$tag}{'tag_E_value'};
			}
		}
	}
=cut
	
##################################### Version 1.09###################################################
# output the best peptide using the combination of multiple tags
	my $peptide_ref;
	my $peptide_ref_name;
	my $i=0;
	my %PeptideEvaluehash;
	
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{
		foreach my $tag (keys %{$sort_results{$pvalue}})
		{
			foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}{$tag}})
			{
				$i++;
				my $j=1;
				foreach my $protein (keys %{$sort_results{$pvalue}{$tag}{$pep}})
				{
					if($j==1)
					{
						my ($M_H, $lSideMass, $rSideMass, $Peptidematchedratio, $TagSeq, $TagEvalue, $tagsidemass, $ions ,$Reference, $Peptide) = (split/\s+/,$sort_results{$pvalue}{$tag}{$pep}{$protein})[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10];

						my ($matched,$total)=split(/\//,$ions);
############## speed #########
						if(!defined($PeptideEvaluehash{$total}{$matched}))
						{
							$PeptideEvaluehash{$total}{$matched} = $self->get_peptide_pvalue($matched,$total);
							$PeptideEvalue = $PeptideEvaluehash{$total}{$matched};
						}
						else
						{
							$PeptideEvalue = $PeptideEvaluehash{$total}{$matched};						
						}
						$sort_results{$pvalue}{$tag}{$pep}{$protein} = sprintf("%-0.3f    %-0.4f    %-0.4f    %-0.2f     %-10s   %-0.2f    %-0.4f  %-10s  %-30s    %-20s", $M_H, $lSideMass, $rSideMass, $PeptideEvalue, $TagSeq, $TagEvalue, $tagsidemass, $ions ,$Reference, $Peptide);
						push @{$peptide_ref->{$pep}->{'info'}}, 
						[$j, $M_H, $lSideMass, $rSideMass, $PeptideEvalue, $TagSeq, $TagEvalue];
						$peptide_ref->{$pep}->{'ref'}->{$protein}++;
						$peptide_ref_name->{$pep}->{'ref'} = $protein;
						$j++;
					}
					else
					{
						$peptide_ref->{$pep}->{'ref'}->{$protein}++;
						$j++;
					}
				}
			
			}
		}
		last if ($i>=$parameter->{'number_of_detailed_result'});			
	}
######################### Output the selected candidate peptide ####################

	my $WeightedEvalue = $self->calculate_weighted_score($peptide_ref);
	print OUTPUT "[Selected identified peptide(s)]\n";
    print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  WeightedEvalue  TagsNum      Reference                       Peptide    \n";
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
	my $order = 0;
    for my $Peptide (sort {
            $WeightedEvalue->{$b}->{'Weighted_Evalue'}
            <=>
            $WeightedEvalue->{$a}->{'Weighted_Evalue'} 
			||
            $peptide_ref_name->{$b}->{'ref'}
            cmp
            $peptide_ref_name->{$a}->{'ref'} 			
		} keys %$WeightedEvalue) #TODO: any better sort sub?
    {
		my $Tag_info    = $WeightedEvalue->{$Peptide}->{'longest_Tag_info'};
        my $each_info   = join "    ", (@$Tag_info)[1, 2, 3, 4, 5];
        my $Evalue      = $WeightedEvalue->{$Peptide}->{'Weighted_Evalue'};
        my $Num_of_Tags = $WeightedEvalue->{$Peptide}->{'Num_of_Tags'};
        my $refs        = $peptide_ref->{$Peptide}->{'ref'};
        my @References  = sort keys %$refs;
        my $first_ref   = shift @References;
		$order++;

		print OUTPUT $order,"      ";
		print OUTPUT $each_info,"           ";
		printf OUTPUT ("%-0.2f           %-2d       %-30s  %-30s\n", $Evalue,$Num_of_Tags,$first_ref,$Peptide);
#        print OUTPUT "$order  $each_info  $Evalue  $Num_of_Tags  $first_ref   $Peptide\n";
        print OUTPUT "                                                                                            $_\n" for @References;
		last if ($order >= $parameter->{'number_of_selected_result'});		
    }	
	
	print OUTPUT "\n\n";
	print OUTPUT "[All identified peptide(s)]\n";	
	print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  TagEvalue  TagSideMass  Ions   Reference                 Peptide   \n";  
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";



######################### Output in detail #####################################	
	my $i=0;
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{
		foreach my $tag (keys %{$sort_results{$pvalue}})
		{
			foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}{$tag}})
			{
				$i++;
				my $j=1;
				foreach my $protein (keys %{$sort_results{$pvalue}{$tag}{$pep}})
				{
					if($j==1)
					{


						printf OUTPUT ("%2d    %10s\n",$i,$sort_results{$pvalue}{$tag}{$pep}{$protein});
						$j++;
					}
					else
					{
						print OUTPUT "                                                                                                 ",$protein,"\n";
						$j++;
					}
				}
				last if ($i>=$parameter->{'number_of_detailed_result'});
			}
		}
	}
	close(OUTPUT);
}

sub WriteResults4NoTags
{
	my ($self,$outfile,$matches,$tag_rank_num) = @_;
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";
	my $parameter = $self->get_parameter();
	my @mod_symbol = ("@","#","%","^","&","*","?","~","!","(",")","{","}","[","]",":",";","'","<",">");	
	
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";


	print OUTPUT "\nJUMP version 2.01 (c) 2012\n\n";
	print OUTPUT "St Jude Children's Research Hospital, X Wang/J Peng\n";
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	printf OUTPUT "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon,$mday,$hour,$min,$sec;	
	my $dtafile = $self->get_scanNum();
	$dtafile =~ s/tag/dta/;
	print OUTPUT "DTA file = ", $dtafile,"\n";
	
	print OUTPUT "Database = ", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass = ",$self->get_precMass(),"   ", "Percentage of precursor peak intensity = ", sprintf("%0.2f",$prec_peak_int*100),"\%\n";
	print OUTPUT "Precursor matched peptides = ", $self->{'matched_prec_peptide_num'},"\, ","Peptides mass tolerance =", $self->get_mass_tolerance(), "  Peptides mass tolerance units = $parameter->{'peptide_tolerance_units'}\n";
	print OUTPUT "Tag number = ",$self->get_tag_number(),"\n";
	
	print OUTPUT "ion series ABCDVWXYZ: ", $parameter->{'ion_series'},"\n";
###### output modification information ###############################
	my $i=0;
	foreach (keys %$parameter)
	{
		if($_ =~/add_/ and $parameter->{$_} != 0)
		{
			print OUTPUT $_,"=",$parameter->{$_},"  ";
		}

		if($_ =~/dynamic_/ and $parameter->{$_} != 0)
		{
			print OUTPUT $_,$mod_symbol[$i],"=",$parameter->{$_},"  ";
			$i++;
		}
	}
	
	print OUTPUT "\n\n";	

	
	my $i=0;
	my %sort_results=();
	my %comb_results=();
	my %peptide_pvalue=();
	my %PeptideEvaluehash=();
	
	foreach my $pep (keys %{$matches})
	{	
		
		my $j=1;
		foreach my $mod (keys %{$matches->{$pep}})
		{	

			my $peptideseq = $matches->{$pep}->{$mod}->{'matchPeptide'};

			if(!$peptide_pvalue{$pep})
			{	

				my $matched = $self->get_peptide_matched($peptideseq,$mod);
				next if ($matched==0);
				my $total = $self->{'Theor_ion_num'};
				my $pep_pvalue = 0;
				if(!defined($PeptideEvaluehash{$total}{$matched}))
				{				
					my $pep_pvalue = $self->get_peptide_pvalue($self->{'Matched_ion_num'} , $self->{'Theor_ion_num'});
					$PeptideEvaluehash{$total}{$matched} = $pep_pvalue;
				}
				else
				{
					$pep_pvalue=$PeptideEvaluehash{$total}{$matched};
				}
				$peptide_pvalue{$pep} = $pep_pvalue;
	#			print $pep,"\t",$mod,"\t",$matched,"\t",$total,"\t",$self->{'Matched_ion_num'},"\t",$self->{'Theor_ion_num'},"\t",$pep_pvalue,"\t",$peptide_pvalue{$pep},"\n";
			}	
			my $ions = $self->{'Matched_ion_num'} . "/$self->{'Theor_ion_num'}";
			
			my $results = sprintf("%-0.3f    %-0.4f    %-0.4f       %-0.2f         %-6s     %-0.2f      %-0.4f       %-6s     %-30s    %-30s", 
			$matches->{$pep}->{$mod}->{'theoretical_MH'}/1000,$matches->{$pep}->{$mod}->{'lthreot_Mass'},$matches->{$pep}->{$mod}->{'rthreot_Mass'},
			$peptide_pvalue{$pep},$matches->{$pep}->{$mod}->{'tagSeq'},$matches->{$pep}->{$mod}->{$parameter->{'tag_select_method'}},$matches->{$pep}->{$mod}->{'sideMass'},$ions,$matches->{$pep}->{$mod}->{'proteinid'},
				$matches->{$pep}->{$mod}->{'matchPeptide_orig'});
			$sort_results{$peptide_pvalue{$pep}}{$tag}{$matches->{$pep}->{$mod}->{'matchPeptide_orig'}}{$matches->{$pep}->{$mod}->{'proteinid'}}=$results;
			
		}
	}


	print OUTPUT "[Selected identified peptide(s)]\n";
    print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  WeightedEvalue  TagsNum      Reference                       Peptide    \n";
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
	my $order = 0;
	my $i=0;
	
	LOOP:
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{
		foreach my $tag (keys %{$sort_results{$pvalue}})
		{
			foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}{$tag}})
			{
				$i++;
				my $j=1;
				foreach my $protein (keys %{$sort_results{$pvalue}{$tag}{$pep}})
				{
					if($j==1)
					{
						my @result_array = split(/\s+/,$sort_results{$pvalue}{$tag}{$pep}{$protein});
						printf OUTPUT ("%2d    %0.3f     %0.4f    %0.4f      %.2f   %10s       %.2f             %1d   %20s  %30s \n",$i,$result_array[0],$result_array[1],$result_array[2],$result_array[3],$result_array[4],$result_array[3],0,$protein,$pep);
				#		printf OUTPUT ("%2d    %10s\n",$i,$sort_results{$pvalue}{$tag}{$pep}{$protein});
						$j++;
					}
					else
					{
						print OUTPUT "                                                                                      ",$protein,"\n";
						$j++;
					}
				}
				last LOOP if ($i>=$parameter->{'number_of_selected_result'});
			}
		}
	}

	print OUTPUT "\n\n";

	
	print OUTPUT "[All identified peptide(s)]\n";	
	print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  TagEvalue  TagSideMass     Ions      Reference                 Peptide   \n";  
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";



######################### Output in detail #####################################	
	my $i=0;
	foreach my $pvalue (reverse sort { $a <=> $b } keys %sort_results)
	{
		foreach my $tag (keys %{$sort_results{$pvalue}})
		{
			foreach my $pep (sort {$a<=>$b} keys %{$sort_results{$pvalue}{$tag}})
			{
				$i++;
				my $j=1;
				foreach my $protein (keys %{$sort_results{$pvalue}{$tag}{$pep}})
				{
					if($j==1)
					{
						printf OUTPUT ("%2d    %10s\n",$i,$sort_results{$pvalue}{$tag}{$pep}{$protein});
						$j++;
					}
					else
					{
						print OUTPUT "                                                                                                  ",$protein,"\n";
						$j++;
					}
				}
				last if ($i>=$parameter->{'number_of_detailed_result'});
			}
		}
	}

	close(OUTPUT);
}


sub calculate_weighted_score # private method
{   
	my ($self, $peptide_ref) = @_;
	my $parameter = $self->get_parameter();	
	my $coeff = $parameter->{'tag_coeff_evalue'};
	my $method=2;
    my $WE;
    ((caller)[0] eq ref $self) or die "Private method called.\n";
    for my $Peptide (keys %$peptide_ref)
    {   
		my $info        = $peptide_ref->{$Peptide}->{'info'};
        #my $Reference   = $peptide_ref->{'ref'};
        #--------
        $WE->{$Peptide}->{'Num_of_Tags'} = scalar @$info;
        my $PeptideEvalue       = $info->[0]->[4]; # all are equal for one peptide
        my @sorted_lines        = sort { length $b->[5] <=> length $a->[5] } @$info;
        my $longest_tag_line    = shift @sorted_lines;
        $WE->{$Peptide}->{'longest_Tag_info'} = $longest_tag_line;
        my ($longest_Tag_seq, $longest_Tag_Evalue) = (@$longest_tag_line)[5, 6];
######## Get the highest score  
		@sorted_lines        = sort { $b->[6] <=> $a->[6] } @$info;
        $highestscore_tag_line    = shift @sorted_lines;
        my ($highest_Tag_seq,$highest_Tag_Evalue) = (@$highestscore_tag_line)[5, 6];
		print $longest_Tag_seq,"\t",$highest_Tag_Evalue,"\n";
#####################################################		
		
        $WE->{$Peptide}->{'Weighted_Evalue'} = $PeptideEvalue + $highest_Tag_Evalue;
        my @TagSeq_Evalue       = map { [$_->[5], $_->[6]] } @sorted_lines;
		#print $PeptideEvalue,"\t",$longest_Tag_Evalue,"\t";
        for my $each_pair (@TagSeq_Evalue)
        {   
			my ($TagSeq, $TagEvalue) = @$each_pair;
        #    $WE->{$Peptide}->{'Weighted_Evalue'} += ($TagSeq =~ /$longest_Tag_seq/) ?  $TagEvalue : ($coeff  * $TagEvalue);
			if($method==1)
			{
				next if ($TagEvalue == $highest_Tag_Evalue);
				next if (length($highest_Tag_seq)>length($TagSeq) and  $highest_Tag_seq=~/$TagSeq/);
				my @substr = $self->lc_substr($highest_Tag_seq,$TagSeq);
				if(scalar(@substr)==0)
				{
					$WE->{$Peptide}->{'Weighted_Evalue'} += $TagEvalue;
				}
				else
				{			
					$WE->{$Peptide}->{'Weighted_Evalue'} += ($TagEvalue == $highest_Tag_Evalue) ?  0 : ($coeff  * $TagEvalue);
				}
			}
			elsif($method==2)
			{
########## compare the other tag sequence with tag with highest score, if the latter sequence already covers the tag sequences, then ignore
####################################################################### if the new tag is longer than the tag with highest score, then calculate a coefficient for the tag score 

				next if ($TagEvalue == $highest_Tag_Evalue);
				next if (length($highest_Tag_seq)>length($TagSeq) and  $highest_Tag_seq=~/$TagSeq/);
				
				my @substr = $self->lc_substr($highest_Tag_seq,$TagSeq);
				if(scalar(@substr)==0)
				{
					$WE->{$Peptide}->{'Weighted_Evalue'} += $TagEvalue;
				}
				else
				{
					$WE->{$Peptide}->{'Weighted_Evalue'} += (length($TagSeq) - scalar(@substr)) / length($highest_Tag_seq) * $TagEvalue;
				}
			
			}
		#	print $WE->{$Peptide}->{'Weighted_Evalue'},"\t",$coeff,"\t", $TagEvalue,"\n";
        }
    }
    return $WE;
}

sub lc_substr {
  my ($self, $str1, $str2) = @_; 
  my $l_length = 0; # length of longest common substring
  my $len1 = length $str1; 
  my $len2 = length $str2; 
  my @char1 = (undef, split(//, $str1)); # $str1 as array of chars, indexed from 1
  my @char2 = (undef, split(//, $str2)); # $str2 as array of chars, indexed from 1
  my @lc_suffix; # "longest common suffix" table
  my @substrings; # list of common substrings of length $l_length
 
  for my $n1 ( 1 .. $len1 ) { 
    for my $n2 ( 1 .. $len2 ) { 
      if ($char1[$n1] eq $char2[$n2]) {
        # We have found a matching character. Is this the first matching character, or a
        # continuation of previous matching characters? If the former, then the length of
        # the previous matching portion is undefined; set to zero.
        $lc_suffix[$n1-1][$n2-1] ||= 0;
        # In either case, declare the match to be one character longer than the match of
        # characters preceding this character.
        $lc_suffix[$n1][$n2] = $lc_suffix[$n1-1][$n2-1] + 1;
        # If the resulting substring is longer than our previously recorded max length ...
        if ($lc_suffix[$n1][$n2] > $l_length) {
          # ... we record its length as our new max length ...
          $l_length = $lc_suffix[$n1][$n2];
          # ... and clear our result list of shorter substrings.
          @substrings = ();
        }
        # If this substring is equal to our longest ...
        if ($lc_suffix[$n1][$n2] == $l_length) {
          # ... add it to our list of solutions.
          push @substrings, substr($str1, ($n1-$l_length), $l_length);
        }
      }
    }
  }   
 
  return @substrings;
}



sub exportMatches
{
	my ($self, $f, $matches) = @_;
#### changed by xusheng ########

	open OUT, '>'.$f or die;
	foreach my $tid (sort keys %$matches)
	{
		foreach my $pep (keys %{$matches->{$tid}})
		{

			next if ($matches->{$tid}->{$pep}->{'tagSeqMatch'}==0);
			print OUT $matches->{$tid}->{$pep}->{'scanNum'},"\t",$matches->{$tid}->{$pep}->{'precMass'},"\t",$matches->{$tid}->{$pep}->{'threotical_MH'},"\t",$matches->{$tid}->{$pep}->{'tagSeq'},"\t",$matches->{$tid}->{$pep}->{'matchPeptide'}->{'seq'},"\t",$pep,"\t"; #$matches->{$tid}->{$pep}->{'matchPeptide'}->{'desc'},"\t";
			if($matches->{$tid}->{$pep}->{'tagSeqMatch'}==1)
			{
				print OUT "Tag sequence matched\t";
			}
			elsif($matches->{$tid}->{$pep}->{'tagSeqMatch'}>1)
			{
				print OUT "Partial tag matched\t";
			}
			else
			{
				print OUT "Tag sequence unmatched\t";
			}
			print OUT $matches->{$tid}->{$pep}->{'sideMass'},"\t";
			if($matches->{$tid}->{$pep}->{'lthreot_Mass'})
			{
				print OUT $matches->{$tid}->{$pep}->{'lthreot_Mass'},"\t";
			}
			else
			{
				print OUT "\t";
			}
				
			if($matches->{$tid}->{$pep}->{'lsideMassMatch'}>0 or $matches->{$tid}->{$pep}->{'rsideMassMatch'}>0)
			{

				print OUT "lSide mass matched\t";
				
			}
			else
			{
				print OUT "lSide mass unmatched\t";
			}
			
			if($matches->{$tid}->{$pep}->{'rthreot_Mass'})
			{
				print OUT $matches->{$tid}->{$pep}->{'rthreot_Mass'},"\t";
			}
			else
			{
				print OUT "\t";
			}
			if($matches->{$tid}->{$pep}->{'rsideMassMatch'}>0)
			{
				print OUT "rSide mass matched\t";
			}
			else
			{
			
				print OUT "rSide mass unmatched\t";
			}
			
			print OUT $matches->{$tid}->{$pep}->{'rankP'},"\t",$matches->{$tid}->{$pep}->{'hyperP'},"\n";
		}
	}
	close OUT;
}	


1;

__END__