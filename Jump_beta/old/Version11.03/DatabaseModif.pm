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

###### The package of task1 has been renamed because there is same fuction in Algorithm::Combinatorics. 
###### USAGE: my %result = changeValues(\%Motif, $comb_num, \%massHash);
#!/usr/bin/perl -I /home/xwang4/scripts
########################### Example #################################
#use strict;
#use warnings;
#use Spiders::DatabaseModif;

#my %modif = ( 'M' => 57, 'T' => 78, 'S' => 79 );
#my %new_modif = ('M'=>'M@','T'=>'T#','S'=>'S$');
#my $comb_num = 6;

#my %pephash = (
#    '1050.25124' => ['PITWGSDVAR'],
#);

#my $cmdatabase = new Spiders::DatabaseModif();
#$cmdatabase->GenerateModifSeq(\%modif,$comb_num,\%pephash,\%new_modif);

#foreach my $key ( keys(%pephash) ) {
#    print "$key=>";
#    foreach ( @{ $pephash{$key} } ) {
#        print "$_,";
#    }
#    print "\n";
#}
#######################################################################

package Spiders::DatabaseModif;
      
use strict;
use warnings;
use Exporter;
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

### This subroutine is used to change the mass value for all possible modifications for database########
### 
sub GenerateModifSeq
{
	my ($self,$modif,$combNUM,$pephash,$new_modif)=@_;
####### generate combination of modified amino acid 
    my %mid_hash =$self->modify_combinations(\%$modif, $combNUM);
############ For each peptide in the peptide hash, change the mass value #######		
    foreach my $pepmass (keys(%$pephash))
    {
        foreach my $peptide (@{$$pephash{$pepmass}})
        {

			my $flag        = 0;
		
			my @modif_array = ();
########## The following codes might not efficient enough for large database 			
			my @peptide_aa = split(//, $peptide);
			my $i=0;
			my %modif_pos_aa;
			foreach my $aa (@peptide_aa)
			{
				$i++;
				if (exists($$modif{$aa}))
				{
					$flag++;
############# save the position, rather than amino acid itself ################					
					$modif_pos_aa{$i}=$aa;			
				}
			}
######### @modif_array is used to save amino acid 			
			@modif_array = map {$_} sort {$a<=>$b} keys %modif_pos_aa;	
####### If there is more than one modification in the peptide ############
			if ($flag != 0)
			{
				for (my $i = 1 ; $i <= $flag ; $i++)
				{

################## do not pass the @array as a parameter (??????????) ###########################					
					my $iter = $self->comb($i,\@modif_array);

					foreach my $c (@$iter)
					{

						my $changed_peptide = $peptide;
						my $i=0;
						foreach my $pos (@$c)
						{
########replace the amino acid with modified amino acid ##############################
								substr($changed_peptide,$pos-1+$i,1,$new_modif->{$modif_pos_aa{$pos}});
								$i++;
						}					
############ @pos_aa is used to save amino acid, rather than position ##################	

						my @pos_aa = map {$modif_pos_aa{$_}} @$c;

						my $join = join("", sort(@pos_aa));

						if (exists($mid_hash{$join}))
						{
				############### mass = modif_mass * 1000		
							my $modif_pepmass = $pepmass + $mid_hash{$join} * 1000;
							
					#		print $modif_pepmass,"\t",$changed_peptide,"\n";
							
							push (@{$pephash->{$modif_pepmass}}, $changed_peptide);
						}							
					}
				}
			}
		}
    }
}

sub modify_combinations
{
	my ($self,$hash,$comb_num) = @_;
    my %new_hash;

    for (my $i = 0 ; $i < $comb_num ; $i++)
    {
        if ($i == 0)
        {
            foreach my $key (keys(%$hash))
            {
                $new_hash{$key} = $$hash{$key};
            }
        }
        else
        {
            foreach my $hash_key (keys(%new_hash))
            {
                foreach my $key (keys(%$hash))
                {
                    my $key_join = "$hash_key" . "$key";
                    my @alpha    = sort(split(//, $key_join));
                    my $key_comb = join('', @alpha);

                    $new_hash{$key_comb} = $$hash{$key}  + $new_hash{$hash_key};

                }
            }
        }
    }
    return (%new_hash);
}

sub comb{
	my ($self, $n, $dataref)=@_;
	my @data=@$dataref;
	my @result;
	return [map {[$_]} @data] if $n==1;
	while(1)
	{
		last if @data<$n;
		my $item=shift @data;
		my $ret=$self->comb($n-1,\@data);
		for(@$ret)
		{
		 unshift @$_,$item;

		 push @result,$_;
		}
	}
	return \@result;
}

1;

