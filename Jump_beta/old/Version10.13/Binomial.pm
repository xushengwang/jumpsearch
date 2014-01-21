#!/usr/bin/perl

######### Bionomial  ##########################################
#                                                             #
#       **************************************************    #  
#       **** Binomial module for Tag		          ****    #     
#       ****					                      ****    #  
#       ****Copyright (C) 20212 - Xusheng Wang	      ****    #     
#       ****all rights reserved.		              ****    #  
#       ****xusheng.wang@stjude.org		              ****    #  
#       ****					                      ****    #  
#       ****					                      ****    #  
#       **************************************************    # 
###############################################################

package Spiders::Binomial;
        
use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

sub new{
	my ($class,%arg)=@_;
    my $self = {
        _n => undef,
		_k => undef,
		_p =>undef,
    };
    bless $self, $class;
	return $self;
}

sub binomial {
	my ($n, $k, $p) = @_;
	if (($p < 0) or (1 < $p)) {
		die "Probabilities must be between 0 and 1"; 
	}
	my $sum=0;
    for ($k = 0; $k <= $n; $k++) {
		return $k == 0 if $p == 0;
		return $k != $n if $p == 1;
		return choose($n, $k) * $p**$k * (1-$p)**($n-$k);
		$sum += $prob;
	}
}

sub choose {
	my ($n, $k) = @_;
	my ($result, $j) = (1, 1);
	return 0 if $k > $n || $k < 0;
	$k = ($n - $k) if ($n - $k) < $k;
	while ( $j <= $k ) 
	{
		$result *= $n--;
		$result /= $j++;
	}
	return $result;
}