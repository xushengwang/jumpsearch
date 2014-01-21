#!/usr/bin/perl

######### Simulation ##########################################
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
# added the normalization function (v1.0.1) on 5/11/2012

package Spiders::Consolidation;
        
use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.01;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
        _dta_file => $arg{'-dta_file'},
		_keepnum =>$arg{'-keepnum'},
    };
    bless $self, $class;
	return $self;
}

sub get_dta
{
	my $self=shift;
	return $self->{_dta_file};
}

sub get_keepnum
{
	my $self = shift;
	return $self->{_keepnum};
}

sub get_strong_peak_window
{
	my $self=shift;
	return $self->{'_strong_peak_window'};
}

sub Consolidation
{
	my $self = shift;

	my $dta = $self->get_dta();
	my $keep = $self->get_keepnum();


	open (DTA, "$dta") || die "can not open the dta file: $dta\n";
	my %dtahash;
	my $line0 = <DTA>;
	my $line1 = <DTA>;
	chomp $line1;
	my ($startmz, $startint) = split(/ /, $line1);
	$dtahash{$startmz}=$startint;
	
	while(<DTA>)
	{
		chomp $_;
#		$_=~s/[[:cntrl:]]+//g;
		my ($mz,$int)=split(/\s+/,$_);
		$dtahash{$mz}=$int;
	}
	close(DTA);
	system(qq(rm $dta));
	open (DTA, ">$dta");

	print DTA "$line0";

	$startmz = sprintf("%f", $startmz);
	my $nextmz = $startmz + 100; 
	my %temp;
	my $num = keys(%dtahash);
	my $i = 0;
	my $j=0;
	my $win_num=0;
# %strong_peak_window used for saving the mass value of the strongest peak within each window	
	my %strong_peak_window;
	foreach my $mz (sort {$a<=>$b} keys %dtahash)
	{
		$j++;

		my $int = $dtahash{$mz};
		if ($mz == 0){ print "Oops....error in converting"; }
	## Get the mass value and save it to the temp hash, but the last number will not put the temp hash	
		if($mz<=$nextmz)
		{
			$temp{$mz}=$dtahash{$mz};
		}
	## goes to next window .. 	

		if($mz>$nextmz)
		{
			if((scalar keys %temp)>1)
			{
				my $k=0;
				my %keepers;
				for my $mz (reverse sort {$temp{$a}<=>$temp{$b}}  keys %temp)
				{
					if($k==0)
					{
						$strong_peak_window{$win_num}=$temp{$mz};
					}
					$keepers{$mz} = $temp{$mz};
					$k++; 
					last if ($k==$keep);
				}
				for my $mz (sort {$a<=>$b} keys %keepers)
				{
					print DTA "$mz $keepers{$mz}\n";
				}
				$win_num++;
			}
			elsif((scalar keys %temp)==1)
			{
				for my $mz (keys %temp)
				{	
					$strong_peak_window{$win_num}=$temp{$mz};	
					print DTA "$mz $temp{$mz}\n";
				}
				$win_num++;
			}			
			%temp=();
			$temp{$mz}=$dtahash{$mz};			
			$nextmz += 100;				
		}
		if(scalar( keys %dtahash) == $j)
		{
			for my $mz (keys %temp)
			{	
				$strong_peak_window{$win_num}=$temp{$mz};	
				print DTA "$mz $temp{$mz}\n";
			}			
		}		
	}

	close(DTA);
	$self->{'_strong_peak_window'} = \%strong_peak_window;	
}

sub Normalize_Consolidation
{
	my ($self)=@_;
	my $strong_peak_window = $self->get_strong_peak_window();
	my $dtafile = $self->get_dta();
	my $dtahash = $self->read_dta_file($dtafile);
	my $keep = $self->get_keepnum();
		
	my $top_peak = $self->get_top_peak($strong_peak_window);
	
	my $i=0;
	foreach my $mz (sort {$a<=>$b} keys %{$dtahash->{'ms2'}})
	{

		my $window_num = int($i/$keep);
		$i++;
		my $strong_peak_within_window = $strong_peak_window->{$window_num};
		next if (!defined($strong_peak_within_window));

		$dtahash->{'ms2'}->{$mz} = $dtahash->{'ms2'}->{$mz}*($top_peak/$strong_peak_within_window)*(1-0.01*$window_num);

	}
	$self->write_dta_file($dtafile,$dtahash);
}

sub get_top_peak
{
	my ($self,$hash)=@_;
	my $max = 0;
	foreach my $value (keys %$hash ) {
		$max = $hash->{$value} if ($hash->{$value} > $max);
	}
	return $max;
}

sub read_dta_file
{
	my ($self,$dtafile) = @_;
	
	open (DTA, "$dtafile") || die "can not open the dta file: $dtafile\n";
	my %dtahash;
	my $line0 = <DTA>;
	my @data=split(/\s+/,$line0);

	$dtahash{'prec_mz'} = $data[0];
	$dtahash{'prec_charge'} = $data[1];

	while(<DTA>)
	{
		chomp $_;
		my ($mz,$int)=split(/\s+/,$_);
		$dtahash{'ms2'}{$mz}=$int;
	}
	close(DTA);
	return \%dtahash;
}

sub write_dta_file
{
	my ($self,$dtafile,$dtahash) = @_;

	if(-e $dtafile)
	{
		system(qq(rm $dtafile));
	}
	open (DTA, ">$dtafile") || die "can not open the dta file: $dtafile\n";
	print DTA $dtahash->{'prec_mz'}," ",$dtahash->{'prec_charge'},"\n";

	foreach my $mz (sort {$a<=>$b} keys %{$dtahash->{'ms2'}})
	{
		my $intensity = sprintf("%.6f",$dtahash->{'ms2'}->{$mz});
		$mz = sprintf("%.6f",$mz);
		print DTA $mz," ",$intensity,"\n";
	}
	close(DTA);
}

1;