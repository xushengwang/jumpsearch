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
# this value is used to remove the precursor ion
# 5/10/2012 
# add the hypergeometric p-value function in this module


package Spiders::Tag;

use strict;
use warnings;
use Spiders::MassUtils;
use Spiders::WilcoxonRankSum;
use Spiders::Hypergeometric;

use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.02;
##### version 1.0.2 #############
# allow specifying the tag_select_method 
###

@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
		_C_value =>undef,
		_H_value=>undef,
		_parameter=>undef,
    };
    bless $self, $class;
	return $self;
}

sub set_C_value
{
	my ($self,$c_value)=@_;
	$self->{_C_value}=$c_value;	
}

sub get_C_value
{
	my $self=shift;
	if(!defined($self->{_C_value}))
	{
		$self->{_C_value}=1.00335;
	}
	return $self->{_C_value};
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


sub set_dta
{
	my ($self,$dta)=@_;
	$self->{'_dta_file'} = $dta;
	return $self->{'_dta_file'};
}


sub get_dta
{
	my $self=shift;
	return $self->{'_dta_file'};
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

sub set_keepnum
{
	my $self = shift;
	return $self->{_keepnum};
}

sub get_keepnum
{
	my $self = shift;
	return $self->{_keepnum};
}



sub get_precursor_mz
{
	my ($self)=@_;
	return $self->{'prec_mz'};
}

# get the msms value from dta file after deisotope and consolidation
sub get_msms_dta
{
	my $self=shift;
	my $H = $self->get_H_value();
	my $msms_hash;
	my $dta = $self->get_dta();
	my $scan;
	if($dta =~ /(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/)
	{
		$scan = $3;
	}
	open(DTA,"$dta") || die "can not open the dta file: $dta\n";
	my $header = <DTA>;
	my ($mz,$charge) = split(/\s+/,$header);
	my $one_charge_mz = ($mz-$H)/$charge+$H;
	$one_charge_mz = sprintf("%.6f", $one_charge_mz);
	$msms_hash->{$scan}->{"prec_mz"} = $one_charge_mz;
	my @mz_array=();
	my @int_array=();
	while(<DTA>)
	{
		my ($mz,$int) = split(/\s+/,$_);
		push(@mz_array,$mz);
		push(@int_array,$int);
	}
	$msms_hash->{$scan}->{'msms_mz'}=\@mz_array;
	$msms_hash->{$scan}->{'msms_int'}=\@int_array;
	close(DTA);
	$self->{'prec_mz'} = $mz;
	return $msms_hash;
}


## get the peaks whose mz are greater than that precursor peaks
sub get_peaks4tags
{
	my ($self,$scan,$msms_hash)=@_;

	#my $init_mz = $prec_mz;
	my $init_mz = 0;
	my $tag_peaks;
	my $mz_array_ref = $msms_hash->{$scan}->{'msms_mz'};
	my $int_array_ref = $msms_hash->{$scan}->{'msms_int'};
	my @mz_array = @$mz_array_ref;
	my @int_array =@$int_array_ref;
	
	for(my $i=0;$i<$#mz_array;$i++)
	{
		if($mz_array[$i]>$init_mz)
		{
			$tag_peaks->{$mz_array[$i]} = $int_array[$i];
		}
	}
	return $tag_peaks;
}

sub derive_tag
{
	my ($self,$dtafile,$msms_hash) = @_;
	my $prec_mz = $self->get_precursor_mz();
	
	my $scan;
	if($dtafile =~ /(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/)
	{
		$scan = $3;
	}
	my $parameter = $self->get_parameter();
	my $tolerance = $parameter->{'tag_tolerance'};
	
	my $mutils = new Spiders::MassUtils();
	$mutils->set_parameter($parameter);

	my $AA_Mass = $mutils->get_AA_mass();

	my $mass_range = $mutils->get_mass_range($AA_Mass);
######### Version 1.03 correct this into 19 
#	my $amino_acid_num = $mutils->get_amino_acid_number($AA_Mass);
	my $amino_acid_num = 19;	
#######################################################	
	
	my $min_mass = $mutils->get_min_mass($AA_Mass);
	my $max_mass = $mutils->get_max_mass($AA_Mass);

	my $tag_peaks = $self->get_peaks4tags($scan,$msms_hash);
	return if (!defined($tag_peaks));

	my $sum_selected_intensity = $self->sum_intensity($tag_peaks);
	my $peaks_rank = $self->rank_intensity($tag_peaks);
	my %tag_AA_hash;

	foreach my $mz (sort {$a<=>$b} keys %$tag_peaks ) 
	{
			
# $tag_AA_hash is used to store all mz and int for those peaks that can be derived tag
		$self->calculate_tag($AA_Mass,$mz,$tag_peaks,\%tag_AA_hash,$min_mass,$max_mass);
		
	}
	
	my %tag_hyper_pvalue_hash;
	$self->calculate_hyper_pvalue($tag_peaks,\%tag_AA_hash,\%tag_hyper_pvalue_hash,$mass_range,$amino_acid_num,$min_mass,$max_mass);
	
	my @tag_array;
	foreach my $mz_pos1 (sort {$a<=>$b} keys %tag_AA_hash)
	{
#		print $mz_pos1,"bb\n"; 
			
		my $tag="";
		my @mz;
		my $position=0;
		
		my %mass_position;
		$self->iterate($mz_pos1,\%tag_AA_hash,\$tag,$position,\%mass_position,\@tag_array,\@mz,$tag_peaks,\%tag_hyper_pvalue_hash);				
	}

	if((scalar @tag_array)>0)
	{
	
########## if the same tag derives from both b+ and y+ ion series, then the E value of tag is the sum of the tag 
=head
		my %temp_hash=();
		foreach (@tag_array)
		{
			if(defined($temp_hash{$_->{'tag'}}))
			{
				if($temp_hash{$_->{'tag'}}>($prec_mz-$_->{'tag_mz'}-$tolerance) and $temp_hash{$_->{'tag'}}<($prec_mz-$_->{'tag_mz'}+$tolerance))
				{
					$temp_hash{$_->{'p_'}}
				}
			}
			else
			{
				$temp_hash{$_->{'tag'}}=$_->{'tag_mz'};
			}
		}
=cut
#################################################################################################################	
	
		my @sorted_tag=();
		if($parameter->{'tag_select_method'} eq "rank_p")
		{
	#sort by rank sum p value
			
			@sorted_tag = sort {$a->{'p_value'}<=>$b->{'p_value'}} @tag_array;
		
		}
		elsif($parameter->{'tag_select_method'} eq "hyper_p")
		{
	#sort by hypergeometric p value	
			@sorted_tag = sort {$a->{'hyper_p'}<=>$b->{'hyper_p'}} @tag_array;	
		}
		elsif($parameter->{'tag_select_method'} ne ["rank_p"|"hyper_p"])
		{
			print "please specify a right parameter for tag_select_method!\n";
		}

#		foreach (@sorted_tag)
#		{
#			print $_->{'tag'},"\t",$_->{'p_value'},"\n";
#		}				
#		$self->select_top_rank_pvalue(\@tag_array);
#		$self->select_top_hyper_pvalue(\@tag_array);
#		return 1;
		return \@sorted_tag;
		
	}
	else
	{
		return 0;
	}
}


sub construct_search_tag
{
	my ($self, $out_filename, $tag) = @_;
	open(OUTPUT,">>$out_filename") || die "can not open the tag output file: $out_filename!";

	my $rank_tag = reverse $tag->{'tag'};
#	my $rank_pvalue = $tag->{'p_value'}>0 ? -10*log($tag->{'p_value'}) :  0; 
	
#	my $rank_hyper_pvalue = $tag->{'hyper_p'}>0 ? -10*log($tag->{'hyper_p'}) :  0;
	my $rank_pvalue = -log($tag->{'p_value'})/log(10);

########This value not used any more ##############
####### it changes back because of low resolution
	my $rank_hyper_pvalue=0;
#	if ($tag->{'hyper_p'}<0.0000000001)
#	{
#		$rank_hyper_pvalue = 10;
#	}
#	else
#	{
		my $hyperp = $tag->{'hyper_p'} =sprintf("%.20f",$tag->{'hyper_p'});
		if($hyperp == 0)
		{
			$rank_hyper_pvalue = 20;
		}
		else
		{
			$rank_hyper_pvalue = -log($hyperp)/log(10);
		}
#	}
######################################	
	my $rank_tag_mz = $tag->{'tag_mz'};

	my $prec_mz = $self->get_precursor_mz();

	if((defined $rank_tag_mz) and (scalar @$rank_tag_mz)>1)
	{
		my @sort_rank_tag_mz = sort {$a <=> $b} (@$rank_tag_mz);

		my $rank_left_tag_mz = $sort_rank_tag_mz[-1];

		my $temp_rank_mz = shift @sort_rank_tag_mz;
		foreach my $mz (@sort_rank_tag_mz)
		{
			$temp_rank_mz .= ";$mz"; 
		}

		my $search_tag = {'scanNum'=>$out_filename,
				 'precMass'=>$prec_mz,
				 'tagSeq'=>$rank_tag,
				 'sideMass'=>$rank_left_tag_mz,
				 'rankP'=>$rank_pvalue,
				 'hyperP'=>$rank_hyper_pvalue
				};
	#	$self->{'search_tag'} = $search_tag;
		print OUTPUT $prec_mz,"\t",$rank_tag,"\t",$rank_left_tag_mz,"\t",$rank_pvalue,"\t",$rank_hyper_pvalue,"\t",$temp_rank_mz,"\n";
		return $search_tag;
	}
	else
	{
		return 0;
	}
	close(OUTPUT);
}

	
	
sub select_top_rank_pvalue
{
	my ($self,$tag_array)=@_;
	my $min_rank_pvalue=1;
	my $select_tag;
	my $select_tag_mz;
	my $select_rank_pvalue;
	my $select_hyper_pvalue;
	
	foreach my $tag (@$tag_array)
	{
		if($tag->{'p_value'}<$min_rank_pvalue)
		{
			my $tag_string = reverse $tag->{'tag'};
			$select_tag = $tag_string;
			$select_tag_mz = $tag->{'tag_mz'};

			$select_rank_pvalue = $tag->{'p_value'}>0 ? -log($tag->{'p_value'})/log(10) :  0;
			$select_hyper_pvalue = $tag->{'hyper_p'}>0 ? -log($tag->{'hyper_p'})/log(10) :  0;

			$min_rank_pvalue = $tag->{'p_value'};
#			print $tag_string,"\t",$tag->{'p_value'},"\t",$tag->{'hyper_p'},"\n";
		}
	}

	
#	print "Tag:","\t",$select_tag,"\t",$select_rank_pvalue,"\t",$select_hyper_pvalue,"\t";
	$self->{'_rank_tag_seq'} = $select_tag;
	$self->{'_rank_pvalue'} = $select_rank_pvalue;
	$self->{'_rank_hyperpvalue'} = $select_hyper_pvalue;
	$self->{'_rank_tag_mz'} = $select_tag_mz;

	
}

sub select_top_hyper_pvalue
{
	my ($self,$tag_array)=@_;
	my $min_rank_pvalue=1;
	my $select_tag;
	my $select_tag_mz;	
	my $select_rank_pvalue;
	my $select_hyper_pvalue;
	
	foreach my $tag (@$tag_array)
	{
		if($tag->{'hyper_p'}<$min_rank_pvalue)
		{
			my $tag_string = reverse $tag->{'tag'};
			$select_tag = $tag_string;
			$select_tag_mz = $tag->{'tag_mz'};
			
			 $select_rank_pvalue = $tag->{'p_value'}>0 ? -log($tag->{'p_value'})/log(10) :  0;
			$select_hyper_pvalue = $tag->{'hyper_p'}>0 ? -log($tag->{'hyper_p'})/log(10) :  0;
			$min_rank_pvalue = $tag->{'hyper_p'};
#			print $tag_string,"\t",$tag->{'p_value'},"\t",$tag->{'hyper_p'},"\n";
		}
	}
#	print $select_tag,"\t",$select_rank_pvalue,"\t",$select_hyper_pvalue,"\n";
	$self->{'_hyper_tag_seq'} = $select_tag;
	$self->{'_hyper_pvalue'} = $select_rank_pvalue;
	$self->{'_hyper_hyperpvalue'} = $select_hyper_pvalue;
	$self->{'_hyper_tag_mz'} = $select_tag_mz;
}


sub OutputTag
{
	my ($self, $out_filename) = @_;
	open(OUTPUT,">$out_filename") || die "can not open the tag output file: $out_filename!";

	my $rank_tag = $self->{'_rank_tag_seq'}; 
	my $rank_pvalue = $self->{'_rank_pvalue'}; 
	my $rank_hyper_pvalue = $self->{'_rank_hyperpvalue'};
	my $rank_tag_mz = $self->{'_rank_tag_mz'};
	my $hyper_tag = $self->{'_hyper_tag_seq'};
	my $hyper_rank_pvalue = $self->{'_hyper_pvalue'}; 
	my $hyper_pvalue = $self->{'_hyper_hyperpvalue'};
	my $hyper_tag_mz = $self->{'_hyper_tag_mz'};
	my $prec_mz = $self->get_precursor_mz();
		
	if((defined $rank_tag_mz) and (scalar @$rank_tag_mz)>1)
	{
		my @sort_rank_tag_mz = sort {$a <=> $b} (@$rank_tag_mz);

		my $rank_left_tag_mz = $sort_rank_tag_mz[-1];

		my $temp_rank_mz = shift @sort_rank_tag_mz;
		foreach my $mz (@sort_rank_tag_mz)
		{
			$temp_rank_mz .= ";$mz"; 
		}
		my @sort_hyper_tag_mz = sort {$a <=> $b} (@$hyper_tag_mz);

		my $hyper_left_tag_mz = $sort_hyper_tag_mz[-1];
		my $temp_hyper_mz = shift @sort_hyper_tag_mz;
		foreach my $mz (@sort_hyper_tag_mz)
		{
			$temp_hyper_mz .= ";$mz"; 
		}

		my $search_tag = {'scanNum'=>$out_filename,
				 'precMass'=>$prec_mz,
				 'tagSeq'=>$rank_tag,
				 'sideMass'=>$rank_left_tag_mz,
				 'rankP'=>$rank_pvalue,
				 'hyperP'=>$rank_hyper_pvalue
				};

	#	$self->{'search_tag'} = $search_tag;
		print OUTPUT $prec_mz,"\t",$rank_tag,"\t",$rank_left_tag_mz,"\t",$rank_pvalue,"\t",$rank_hyper_pvalue,"\t",$temp_rank_mz,"\n";
		return $search_tag;
	}
	else
	{
		return 0;
	}
	#$temp_rank_mz,"\t",$hyper_tag,"\t",$hyper_rank_pvalue,"\t",$hyper_pvalue,"\t",$hyper_left_tag_mz,"\t",$temp_hyper_mz,"\n";
}


## This subroutine is difficult to write
# Here is example 
#first pos:955.444885253906      second peak:1069.48815917969
# N
#first pos:1069.48815917969      second peak:1182.56848144531
#NL
#first pos:1182.56848144531      second peak:1311.60998535156
#NLE
#first pos:1311.60998535156      second peak:1412.66735839844
#NLET
#peak: 1182.56848144531 1329.62426757812
#NLE (?)     
# peak: 1182.56848144531 1369.63195800781
#NL (?)

sub iterate
{
	my ($self,$pos1,$tag_AA_hash,$tag,$position,$mass_position,$tag_array,$mz,$tag_peaks,$tag_hyper_pvalue_hash)=@_;
	my $flag = 0;
	foreach my $pos (sort {$a<=>$b} keys %{$tag_AA_hash})
	{
#		next if($pos<=$pos1);
		if (defined($tag_AA_hash->{$pos}->{$pos1}))
		{
		
### to difine the first amino acid .... i.e. N has to control ...		
			if(defined($mass_position->{$pos1}))
			{
				my $seq=substr($$tag,0,$mass_position->{$pos1}-1);
			
				splice(@$mz,2*length($seq));
				$$tag =$seq;
#				$position=0;
#				$$hyper_p=1;
			}		
			$position ++;			
			my $aa = $tag_AA_hash->{$pos}->{$pos1};
			
			my $hyper_pvalue = $tag_hyper_pvalue_hash->{$pos}->{$pos1};

			$$tag = $$tag . $aa;

				
#			$$hyper_p = $$hyper_p * $hyper_pvalue;
#			$$hyper_p = sprintf("%.10f",$$hyper_p);
			

			push(@$mz,$pos1);
			push(@$mz,$pos);
#			print @$mz,"\t",$$tag,"\n";
#			print $$tag,"\n";
			my $hyper_p = 1;
			for(my $i=0;$i<=$#$mz;$i++)
			{
				$hyper_p = $hyper_p * $tag_hyper_pvalue_hash->{$$mz[$i+1]}->{$$mz[$i]};
				$i++;
			}
			my %tag_temp;
			$tag_temp{'tag'}=$$tag;
			$tag_temp{'hyper_p'}=$hyper_p;

			my @temp_mz = @$mz;
			my %temp_hash   = map { $_ => 1 } @temp_mz;
			my @temp_unique = sort keys %temp_hash;			
			$tag_temp{'tag_mz'} = \@temp_unique;

			my $p_value = $self->scoring_tag($mz,$tag_peaks);
			
			$tag_temp{'p_value'} = $p_value;
						
#			print $pos,"\t",$pos1,"\t",$aa,"\t",$hyper_p,"\t",$hyper_pvalue,"\t",$$tag,"\t",$p_value,"\t",$position,"ccc\n";
						
			push (@$tag_array,\%tag_temp);

			$flag = 1;
			$mass_position->{$pos1}=length($$tag);	
			$self->iterate($pos,$tag_AA_hash,$tag,$position,$mass_position,$tag_array,$mz,$tag_peaks,$tag_hyper_pvalue_hash);
		
		}
		else
		{
			if(defined($mass_position->{$pos1}) )
			{

				my $seq=substr($$tag,0, $mass_position->{$pos1}-1);
#				print $$tag,"\t",$seq,"\t",$pos1,"\t",$mass_position->{$pos1},"bbbb\n";				
#				last if (length($seq)<3);
				$$tag =$seq;
	#			splice(@$mz,2*($mass_position->{$pos1}-1));
				$position=length($$tag);
			
			}				
		}		
	}	
	if($flag==0)
	{
		return 0;
	}
}

sub scoring_tag
{
######## scoring using wilcoxon rank sum method
	my ($self,$mz,$tag_peaks)=@_;

	my @sel_array;
	my @unsel_array;
	my $mz_size = scalar(@$mz);
	if($mz_size==2)
	{
		push(@sel_array,$tag_peaks->{$$mz[0]});
		push(@sel_array,$tag_peaks->{$$mz[1]});
		foreach my $unsel_mz (keys %$tag_peaks)
		{

			next if($tag_peaks->{$unsel_mz} eq $tag_peaks->{$$mz[0]});
			next if($tag_peaks->{$unsel_mz} eq $tag_peaks->{$$mz[1]});
			push(@unsel_array,$tag_peaks->{$unsel_mz});
		}
	}
	else
	{
		for(my $i=0;$i<$mz_size;$i++)
		{
			if($i<2)
			{
				push(@sel_array,$tag_peaks->{$$mz[0]});
				push(@sel_array,$tag_peaks->{$$mz[1]});							
			}
			else
			{
				push(@sel_array,$tag_peaks->{$$mz[$i+1]});								
			}
			$i++;
		}
		foreach my $unsel_mz (keys %$tag_peaks)
		{
			foreach my $sel (@sel_array)
			{
				next if($sel eq $tag_peaks->{$unsel_mz});
				push(@unsel_array,$tag_peaks->{$unsel_mz});
			}
		}
	}
	my $wilcox_test = Spiders::WilcoxonRankSum->new();

	$wilcox_test->load_data(\@sel_array, \@unsel_array);
	my $prob = $wilcox_test->probability();
	

	return $prob;
#############			

}

sub calculate_hyper_pvalue
{
	my ($self,$tag_peaks,$tag_AA_hash,$tag_hyper_pvalue_hash,$mass_range,$amino_acid_num,$min_mass,$max_mass) = @_;
########## Version 1.03 for low resolution, the mass tolerance window has to be changed ################### 	
#	my $mass_tolerance = 0.02;
	my $parameter = $self->get_parameter();
	my $mass_tolerance = $parameter->{'tag_tolerance'};
	
	my $mass_poss_bin = int($mass_range/($mass_tolerance*2));
	my $hyper = new Spiders::Hypergeometric();
	
	foreach my $pos (sort {$a<=>$b} keys %{$tag_AA_hash})
	{
######
# 	$min_mass = 57  and $max_mass = 187 if there is no mod
#   $pos = 
		my $number_peak_window = $self->number_peaks($pos,$tag_peaks,$min_mass,$max_mass);			
		my $number_matched_peaks = scalar keys %{$tag_AA_hash->{$pos}};
#		print $mass_poss_bin,"\t",$number_peak_window,"\t",$amino_acid_num,"\t",$number_matched_peaks,"ss\t";
#######
### $mass_poss_bin


		my $hyper_pvalue=$hyper->Hypergeometric($mass_poss_bin,$number_peak_window,$amino_acid_num,$number_matched_peaks);
		
		foreach my $pos1 (sort {$a<=>$b} keys %{$tag_AA_hash->{$pos}})
		{
			$tag_hyper_pvalue_hash->{$pos}->{$pos1} = $hyper_pvalue;
#			print $pos,"\t",$pos1,"\t",$hyper_pvalue,"\n";
		}
	}
		
}

sub number_peaks
{
	my ($self,$pos,$tag_peaks,$min_mass,$max_mass) = @_;
	my $num_peaks=0;
	
	foreach my $mz (keys %$tag_peaks)
	{
## one more unit for counting peaks	
		if($mz<($pos-$min_mass-1) and $mz>($pos-$max_mass+1))
		{
			$num_peaks++;
		}
	}
	return $num_peaks;
}

sub calculate_tag
{
	my ($self,$AA_Mass,$selected_mz,$tag_peaks,$tag_hash,$min_mass,$max_mass)=@_;
	my $tag_mass_error;
	
	my $hyper_geometric = new Spiders::Hypergeometric();
	
	my $parameter = $self->get_parameter();
	my $tolerance = $parameter->{'tag_tolerance'};

	foreach my $mz (keys %$tag_peaks)
	{
		next if($mz eq $selected_mz);
		my $mass_diff = $selected_mz-$mz;
		next if ($mass_diff<0 || abs($mass_diff)>$max_mass || abs($mass_diff)<$min_mass);

		foreach my $aa (keys %$AA_Mass) {
			if((abs($AA_Mass->{$aa} - abs($mass_diff)))<$tolerance)
			{
#				print $aa,"\n";
				
				
				
				$tag_hash->{$selected_mz}->{$mz}=$aa;

				
				$tag_mass_error->{$selected_mz}->{$mz}=abs($AA_Mass->{$aa} - abs($mass_diff));
			}
		}

	}
}


sub sum_intensity
{
	my ($self,$int_hash) =@_;
	my $sum_int=0;
	foreach my $mz (keys %$int_hash ){      
			$sum_int+=$int_hash->{$mz};
	}    
	return $sum_int;
	
}

sub rank_intensity
{
	my ($self,$tag_peaks) =@_;
	my $peaks_rank;
	my $rank=0;
	foreach my $mz (reverse sort {$tag_peaks->{$a}<=>$tag_peaks->{$b}} keys %$tag_peaks ){ 
		$rank++;
		$peaks_rank->{$mz}=$rank;
	}
	return $peaks_rank;
}

1;