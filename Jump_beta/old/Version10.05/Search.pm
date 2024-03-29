package Spiders::Search;

######### Simulation ##########################################
#                                                             #
#       **************************************************    #
#       **** Search database	          ****    #
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
######### version 2.0.1 ############
# change the search module to allow search for low resolution data 
# 1. added the fragmented ion mass tolerance  
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
	my $precmass = $self->get_precMass();
#	my $H = $self->get_H_value();	
	
#	my $max_modif = 0;
######### version 1.1 using index to search the value #####################	
# get precursor mass 
# if there is no modification 
####### ----------------------o-------------------------o------------
####### ---------------------start---------------------end------------
#############################=============selected=======############
	
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
				
				
				while($nomod_pep =~ /[^a-zA-Z]+/g)
				{
					$mod_pos{$-[0]}=1;
					$nomod_pep =~ s/[^a-zA-Z]+//;
				}				
				if((scalar keys %mod_pos)==0)
				{
					my $length=length($peptidehash->{$PeptideID}->{'seq'});
					$nomod_pep=$peptidehash->{$PeptideID}->{'seq'};
					$mod = ":" x $length;
				}
				else
				{
############## if the peptide has modification, get the modification info ##################				
			
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
	#			print $peptidehash->{$PeptideID}->{'seq'},"\t",$nomod_pep,"\t",$mod,"\t";
#				my $peptide_mass = 0;				
				
#				$peptide_mass = $massutil->get_peptide_mass($nomod_pep,$mod,"N");
				
	#			print $prec_mass,"\t",$peptide_mass,"\n";
	
########### In version2.0.1, we merged all individual hashes into one MatchedMass hash 	
				$MatchedMass{$PeptideID}{$mod}{'proteinID'} = $peptidehash->{$PeptideID}->{'proteinID'};
				$MatchedMass{$PeptideID}{$mod}{'seq'} = $peptidehash->{$PeptideID}->{'seq'};
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

#			my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq,$mod,$tags,$error);
			my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq,$mod,$tags,$frag_tolerance);

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

				$matchResults->{$pep}->{$mod}->{'matchPeptide_orig'} = $MatchedMass{$pep}{$mod}{'seq'};

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
							my ($tagSeqMatchResult,$sideMassMatchResult) = $self->matchTagSeq($pepseq,$mod,\%tags_temp,$frag_tolerance);
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
								$matchResults->{$pep}->{$mod}->{'matchPeptide_orig'} = $MatchedMass{$pep}{$mod}{'seq'};							
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
		last if($whole < 0);
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


	
### commented by xusheng; only error without multiple with $tag->{'sideMass'}
	
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

	$result->{'lthreot_Mass'} = $massutil->get_peptide_mass($lSide,$lSide_mod,"C");

	$result->{'rthreot_Mass'} = $massutil->get_peptide_mass($rSide,$rSide_mod,"N");


	
### commented by xusheng; only error without multiple with $tag->{'sideMass'}
	
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
	
	$massutils->getFragmentMasses(pept=>$peptide, modif=>$modif, fragTypes=>['b','y','b++','y++','immo'], spectrum=>\%spectrum);

	my (@aa, @bb);
	
	foreach (%{$spectrum{'mass'}{'term'}}){
		push (@aa, @{$spectrum{'mass'}{'term'}{$_}});
	}
	foreach (%{$spectrum{'mass'}{'intern'}}){
		@bb =  @{$_->{'masses'}};
	}

	my @big = (@aa,@bb);
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
				$matched++;
			}
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

sub get_peptide_pvalue
{
	my ($self,$peptide,$modif)=@_;
	my $theor_mz = $self->generate_peptide_theoritical_mass($peptide,$modif);
	my $exp_mz = $self->get_exp_mz();
	my $matched = $self->compare_theoritical_experiment($exp_mz,$theor_mz);
	
	my $mass_tolerance = $self->get_mass_tolerance();
	
	my $total_location = $#$theor_mz/$mass_tolerance;
	
	my $theor_mass_num = scalar (@$theor_mz);
	my $exp_mass_num = scalar (@$exp_mz);
	
	my $hyper = new Spiders::Hypergeometric();
	my $peptide_pvalue=$hyper->Hypergeometric($total_location,$theor_mass_num,$exp_mass_num,$matched);
	
	my $log_peptide_pvalue = $peptide_pvalue>0 ? sprintf("%.4f",-log($peptide_pvalue)) : 0;

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
	print OUTPUT "Input file = ", $self->get_scanNum(),"\n";
	print OUTPUT "Database =", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass =",$self->get_precMass(),"\n\n";
	close(OUTPUT);
	
}

sub WriteResults
{
	my ($self,$outfile,$matches,$tag_rank_num) = @_;
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";
	my $parameter = $self->get_parameter();
## multiple tag #######	
	open(OUTPUT,">$outfile") || die "can not open the output file: $output!\n";


	print OUTPUT "\nspTag version 2.01 (c) 2012\n\n";
	print OUTPUT "St Jude Children's Research Hospital, X Wang/J Peng\n";
	print OUTPUT "Input file = ", $self->get_scanNum(),"\n";
	print OUTPUT "Database =", $parameter->{'database_name'}, "\n";
	print OUTPUT "Precursor mass =",$self->get_precMass(),"\n\n";
	
#	print OUTPUT "Tag = ", $self->get_tagSeq(),"\n";
#	my $tag_e_value = sprintf("%.2f",$self->get_tagPvalue());
#	print OUTPUT "Tag E_value = ", $tag_e_value,"\n","Tag Rank Number = ",$tag_rank_num,"\n","Tag Side Mass = ",$self->get_sideMass(),"\n"; 
	print OUTPUT "Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  TagEvalue  TagSideMass     Reference                 Peptide   \n";  
	print OUTPUT "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
	my $i=0;
	my %sort_results;
	foreach my $tag (keys %{$matches})
	{
		foreach my $pep (keys %{$matches->{$tag}})
		{	
			foreach my $mod (keys %{$matches->{$tag}->{$pep}})
			{	
	#		print $matches->{$pep}->{'tagSeqMatch'},"\n";
				next if ($matches->{$tag}->{$pep}->{$mod}->{'tagSeqMatch'} < 1);
		#		print $matches->{$pep}->{$mod}->{'lsideMassMatch'},"\t",$matches->{$pep}->{$mod}->{'rsideMassMatch'},"\n";
	#			if ((abs($matches->{$pep}->{$mod}->{'lthreot_Mass'} - $self->get_sideMass()) < 0.02) || (abs($matches->{$pep}->{$mod}->{'rthreot_Mass'} - $self->get_sideMass()) < 0.02))
	#			{
					my $pep_pvalue = $self->get_peptide_pvalue($matches->{$tag}->{$pep}->{$mod}->{'pepseq'},$mod) + $matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}};
	# version 1.02 $matches->{$tag}->{$pep}->{$mod}->{'rank_p'},
	# version 1.03 $matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}}
					my $results = sprintf("%-0.3f    %-0.4f    %-0.4f    %-0.2f     %-10s   %-0.2f    %-0.4f    %-30s    %-30s", 
					$matches->{$tag}->{$pep}->{$mod}->{'theoretical_MH'}/1000,$matches->{$tag}->{$pep}->{$mod}->{'lthreot_Mass'},$matches->{$tag}->{$pep}->{$mod}->{'rthreot_Mass'},
					$pep_pvalue,$matches->{$tag}->{$pep}->{$mod}->{'tagSeq'},$matches->{$tag}->{$pep}->{$mod}->{$parameter->{'tag_select_method'}},$matches->{$tag}->{$pep}->{$mod}->{'sideMass'},$matches->{$tag}->{$pep}->{$mod}->{'proteinid'},
					$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'});
					$sort_results{$pep_pvalue}{$tag}{$matches->{$tag}->{$pep}->{$mod}->{'matchPeptide_orig'}}{$matches->{$tag}->{$pep}->{$mod}->{'proteinid'}}=$results;
					
	#			}
			}
		}
	}	
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
						print OUTPUT "                                                                          ",$protein,"\n";
						$j++;
					}
				}
			}
		}
	}
	close(OUTPUT);
}



sub exportMatches{

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