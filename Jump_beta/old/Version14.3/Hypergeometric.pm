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

package Spiders::Hypergeometric;
        
use strict;
use bignum;

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
		_r =>undef,
    };
    bless $self, $class;
	return $self;
}

sub Hypergeometric {
	my ($self,$n, $k, $r,$x) = @_;
	return $k == 0 if $r == 0;
	return choose($r,$x)*choose($n-$r,$k-$x)/choose($n, $k);
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

1;