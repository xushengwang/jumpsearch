#!/usr/bin/perl

######### Job ##########################################
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
package Spiders::Job;


use strict;
use warnings;
use File::Basename;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.02;
######### Version 1.0.2 #######################
# change the code to allow multiple tags for searching 
#
############################################ 

@ISA	 = qw(Exporter);
@EXPORT      = ();
 

sub new{
    my ($class,%arg) = @_;
    my $self = {
        _dta_path => undef,
        _sim_path  => undef,
    };
    bless $self, $class;
    return $self;
}

sub set_dta_path
{
	my ($self,$dta_path)=@_;
	$self->{_dta_path}=$dta_path;
}

sub get_dta_path
{
	my $self=shift;
	return $self->{_dta_path};
}

sub create_script
{
	my ($self)=shift;
	my $dir = $self->get_dta_path();
	
	open(RUNSHELL,">$dir/runsearch_shell.pl");
print RUNSHELL <<EOF;	
#!/usr/bin/perl  -I /home/xwang4/scripts -I /usr/local/lib/perl5
use Getopt::Long;

use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::Params;
use Spiders::Deisotope;
use Spiders::Consolidation;
use Spiders::Simulation;
use Spiders::ProcessingRAW;
use Spiders::ProcessingMzXML;
use Spiders::Decharge;
use Spiders::Sequest;
use Spiders::Tag;
use Spiders::Error;
use Spiders::Search;
use Spiders::Dta;

my (\$help,\$parameter,\$sim_path);
GetOptions('-help|h'=>\\\$help,
		'-param=s'=>\\\$parameter,
		'-dta_path=s'=>\\\$dta_path,
		);
my \@dtafiles = \@ARGV;		
my \$p = Spiders::Params->new('-path'=>\$parameter);
my \$params=\$p->parse_param();	

my \$mass_tolerance = \$params->{'peptide_tolerance'};
my \$frag_tolerance = \$params->{'frag_mass_tolerance'};
my \$tag_search_method = \$params->{'tag_search_method'};
my \$max_number_tags_for_search = \$params->{'max_number_tags_for_search'};

my \$databasename = basename(\$params->{'database_name'});

foreach my \$dta_file (\@dtafiles)
{	

	next if(\$dta_file eq "\.");
	next if(\$dta_file eq "\.\.");
	next if(\$dta_file !~/\.dta/);

	my \$dta = new Spiders::Dta();
	
#	my \$decharge = new Spiders::Decharge();
#	\$decharge->set_parameter(\$params);
#	\$decharge->set_dta_path(\$dta_path);	
#	my \$dta_hash = \$decharge->create_dtahash(\$dta_file,\\\%msms_hash);
#	my \$dta_file = \$decharge->decharge(\$dta_file, \$dta_hash, \\\%msms_hash, \\\%ms_hash, \\\@mz_array, \\\%realcharge);
	
	my \$deisotope = new Spiders::Deisotope();
	\$deisotope->set_parameter(\$params);
	my \$dtafile = "\$dta_file";
#	my \$dtafile = abs_path(\$dtafile);

	
	\$deisotope->set_dta(\$dtafile);
	\$deisotope->MS2_deisotope();
#	\$deisotope->print_mass_error("\$dtafile.mserr");
	
	my \$consolidation = new Spiders::Consolidation('-dta_file'=>\$dtafile,'-keepnum'=>6);
	\$consolidation->Consolidation();
	\$consolidation->Normalize_Consolidation();
	
	my \$sequest= new Spiders::Sequest();
	\$sequest->Run_Sequest_standalone(\$dtafile);

	my \$tag = new Spiders::Tag();
	\$tag->set_parameter(\$params);
	\$tag->set_dta(\$dtafile);
	my \$msms_hash_dta=\$tag->get_msms_dta(\$dtafile);
	my \$prec_mass = \$tag->get_precursor_mz();
	my \$cand_tag = \$tag->derive_tag(\$dtafile,\$msms_hash_dta);
	my \$tag_rank_num = 0;

    my \$tagfile = \$searchfile = \$dtafile;
    \$tagfile =~ s\/dta\/tag\/;
    \$searchfile =~ s\/dta\/spout\/;	
	
    my \@searchtagarray;	
	foreach \$sel_tag (\@\$cand_tag)
	{	
		\$tag_rank_num++;
		my \$search_tag = \$tag->construct_search_tag(\$tagfile,\$sel_tag);			
        push (\@searchtagarray,\$search_tag);		
	}

	my \$search = new Spiders::Search();
	\$search->set_parameter(\$params);
	
	\$search->set_database("\$dta_path\/.database\/\$databasename");
    \$search->set_precMass(\$prec_mass);

    \$search->set_mass_tolerance(\$mass_tolerance);
    \$search->set_frag_tolerance(\$frag_tolerance);
    my \$searchmassresults = \$search->SearchMass();
    my \$updatedresults;
	if($tag_search_method == 1)
	{
		foreach \$search_tag (\@searchtagarray)
		{
			\$search->set_tag(\$search_tag);
			my \$searchresults = \$search->(\$searchmassresults);
			if(\$searchresults)
			{
				\$updatedresults = \$search->mergeresults(\$updatedresults,\$searchresults,\$search_tag->{'tagSeq'});
				last;
			}		
		}	
	}
	elsif($tag_search_method == 2)
	{
		my $number4search = $max_number_tags_for_search <= $#searchtagarray ? $max_number_tags_for_search : $#searchtagarray;
		
		for(my $i=0;$i<$number4search;$i++)
		{
			\$search->set_tag(\$searchtagarray[$i]);
			my \$searchresults = \$search->(\$searchmassresults);
			if(\$searchresults)
			{
				\$updatedresults = \$search->mergeresults(\$updatedresults,\$searchresults,\$search_tag->{'tagSeq'});
			}		
		}
	}		
	\$dta->set_dta_file(\$dtafile);
	\$dta->process_dtafile();
	my \$exp_mz = \$dta->get_mz_array();
	\$search->set_exp_mz(\$exp_mz);
	\$search->WriteResults(\$searchfile,\$updatedresults,\$tag_rank_num);
}	
EOF
}

sub create_job_files
{
	my ($self,$dta_path)=@_;
	
#	open(DECHARGE,">job.pl");
}

1;
