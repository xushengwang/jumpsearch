#!/usr/bin/perl  -I /home/xwang4/scripts/dev/ 

################################################################
#                                                              #
#       **************************************************     # 
#       ****         JUMP                             ****     # 
#       ****  version 1.16                            ****     # 
#       ****  Xusheng Wang / Junmin Peng              ****     # 
#       ****  Copyright (C) 2012 - 2013               ****     # 
#       ****  All rights reserved.                    ****     # 
#       ****xusheng.wang@stjude.org                   ****     # 
#       ****                                          ****     # 
#       ****                                          ****     # 
#       **************************************************     # 
#                                                              #
################################################################

use Getopt::Long;

use Cwd;
use Cwd 'abs_path';
use Storable;
use File::Basename;
use Spiders::Params;
use Spiders::ProcessingRAW;
use Spiders::ProcessingMzXML;
use Spiders::Decharge;
use Spiders::BuildIndex;
use Spiders::Job;
use Spiders::Error;
use Spiders::Digestion;
use Spiders::MakeDB;
#use Spiders::SummaryResults;
use Spiders::MassAccuracy;
use Spiders::PIP;
use Spiders::SpoutParser;

my $library = "/home/xwang4/scripts/dev/";
my $VERSION = 1.01;

my $progname = $0;
# remove path from our name
$progname =~ s@(.*)/@@i;

my ($help,$parameter,$raw_file);
GetOptions('-help|h'=>\$help,
			'-param=s'=>\$parameter,
			'-rawfile=s'=>\$raw_file,
		);


######## predefine the input #####################
#my $raw_file=$ARGV[0];
usage() if ($help || !defined($parameter) || !defined($raw_file));
check_input($raw_file,\$parameter);
######### programming starting information #############

print <<EOF;
	
################################################################
#                                                              #
#       **************************************************     # 
#       ****                                          ****     # 
#       ****                 JUMP                     ****     # 
#       ****              version 1.16                ****     # 
#       ****        Xusheng Wang / Junmin Peng        ****     # 
#       ****         Copyright (C) 2012 - 2013        ****     # 
#       ****            All rights reserved           ****     # 
#       ****                                          ****     # 
#       **************************************************     # 
#                                                              #
################################################################
EOF

print "\n\n  Initializing JUMP program\n";

###### Get working directory #####################
my $curr_dir = getcwd;

my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();


my $proc_raw = new Spiders::ProcessingRAW();
$proc_raw->set_raw_file($raw_file);
my $dta_path = $proc_raw->get_rawfile_basename();


if(-e($dta_path))
{
		print "  Do you want to remove the old data folder with same name? (y/n): ";
		chomp(my $choice = <STDIN>);
		if($choice eq "yes" || $choice eq "y")
		{	
			system(qq(rm -rf $dta_path >/dev/null 2>&1));
			system(qq(mkdir $dta_path  >/dev/null 2>&1));
		}
}
else
{
			system(qq(mkdir "$curr_dir/$dta_path"  >/dev/null 2>&1));	
}

$p->generateSeqeust($params,'sequest.params');



########################################################################
#my $fasta_file = $params->{'database_name'};
my $databasename = $params->{'database_name'};
my $database_path = dirname($databasename); 
my $database_basename = basename($databasename);

my $num_dynamic_mod = 1;
my $partial_job = 1;
$partial_job++ if ($params->{'digestion'} eq 'partial');

foreach my $key (keys %$params)
{
	if($key=~/dynamic_/)
	{
		$num_dynamic_mod++;
	}
}

if ($params->{search_engine} eq 'SEQUEST')
{
}
elsif(!(-e ($databasename) and $databasename=~/.mdx/))
{
	print "  Creating database\n";
	my $mass_index = $databasename . ".mdx";

	my $peptide_index = $databasename . ".pdx";
	my $protein_index = $databasename . ".prdx";
	my $summary_index = $databasename . ".sdx";
	my $reverse_data = $databasename . ".fasta_reverse";
	if(-e($mass_index))
	{
		print "  Do you want to remove the old database with same name? (y/n): ";
		chomp(my $choice = <STDIN>);
		if($choice eq "yes" || $choice eq "y")
		{
			print "  Removing old database\n";
			system(qq(rm -f $mass_index >/dev/null 2>&1));
			system(qq(rm -f $peptide_index >/dev/null 2>&1));
			system(qq(rm -f $protein_index >/dev/null 2>&1));
			system(qq(rm -f $summary_index >/dev/null 2>&1));
			system(qq(rm -f $reverse_data >/dev/null 2>&1));
		}	
	}
	$params->{'database_name'} = $mass_index;
	################# create database using cluster system ########################
	my $total_protein_number=0;
	open(FASTA, $databasename) || die "can not open the database\n";
	while(<FASTA>)
	{
		$total_protein_number++ if($_=~/^\>/);
	}		
	close(FASTA);
	
#	my $num_protein_per_job = int(4000/($num_dynamic_mod));
	my $num_protein_per_job = int($total_protein_number/(100*$num_dynamic_mod*$partial_job))+1;

	my $protein_num=0;
	my $k=0;
	
	open(FASTA, $databasename) || die "can not open the database\n";
	while(<FASTA>)
	{
		$protein_num++ if($_=~/^\>/);
		if(($protein_num % $num_protein_per_job)==0)
		{
			if($_=~/^\>/)
			{
				$k++;
				print "\r  Generating $k temporary files";
				open(FASTATEMP,">$database_path/temp_${k}_${database_basename}");
			}
			print FASTATEMP "$_";

		}
		else
		{
			print FASTATEMP "$_";		
		}
	}
	print "\r  Generating $k temporary files";	
	close(FASTA);
	print "\n";
	
########### make job file ##############	
	my $job = new Spiders::Job();
	my $abs_parameter = abs_path ($parameter);
	$job->set_library_path($library);
#	my $num_mass_region = 20 * $num_dynamic_mod;
	my $num_mass_region = 20;
	$num_mass_region = $num_mass_region * $partial_job * 2;
	$job->make_createdb_script("$curr_dir/$dta_path",$abs_parameter,$num_mass_region,$num_protein_per_job);

######### submit job file ##############	
	my $job_num=int($protein_num/$num_protein_per_job)+1;	
	for($i=1;$i<=$job_num;$i++)
	{
		open(JOB,">$curr_dir/$dta_path/job_db_$i.sh") || die "can not open the job db files\n";
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q normal\n";
			print JOB "#BSUB -M 20000\n";
			print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
			print JOB "#BSUB -eo $curr_dir/$dta_path/$i.e\n";
			print JOB "#BSUB -oo $curr_dir/$dta_path/$i.o\n";
			print JOB "perl $curr_dir/$dta_path/create_db.pl $database_path/temp_${i}_${database_basename}\n";
			print JOB "rm -f $database_path/temp_${i}_${database_basename}";
			close(JOB);
			system(qq(cd $dta_path && bsub <job_db_$i.sh >/dev/null 2>&1));	
			
		}
		if($params->{'Job_Management_System'} eq 'SGE')
		{
			print JOB "#!/bin/bash\n";
			print JOB "#\$ \-N job_db_$i\n";
			print JOB "#\$ \-e $curr_dir/$dta_path/$i.e\n";
			print JOB "#\$ \-o $curr_dir/$dta_path/$i.o\n";
			print JOB "perl $curr_dir/$dta_path/create_db.pl $database_path/temp_${i}_${database_basename}\n";
			print JOB "rm -rf $database_path/temp_${i}_${database_basename}\n";
			close(JOB);
			system(qq(qsub -cwd -pe mpi 4 -l mem_free=4G,h_vmem=8G "$dta_path/job_db_$i.sh" >/dev/null 2>&1));
			
		}
		close(JOB);
	}
	
	print "  You submit $job_num jobs for creating index files\n";
	
	Check_Job_stat("job_db_",$job_num);

###### Merge all temp_ files into big data file #######################
	
	my $j=0;
	my $l=0;
	my $prev_run = 0;
	print "\n  Sorting indexes\n";

	$job->make_partialidx_script("$curr_dir/$dta_path");
	
	my $prot_job_num = int($num_mass_region/2);

	my $inclusion_decoy = $params->{'inclusion_decoy'};
	
	for($m=0;$m<$num_mass_region;$m++)
	{ 		
		#print $m,"\n";
		my ($mass_hash,$peptide_hash,$protein_hash);	
		#create bash files
		my $parameters = {'range'=> $m,				
						  'GridType'=> $params->{'Job_Management_System'},
				          'JobNum'=>$job_num,
				          'database_path'=> $database_path,
				          'dta_path'=> $dta_path,
				          'mass_index'=>  $mass_index.".$m",
				          'peptide_index'=> $peptide_index.".$m",
				          'protein_index'=> $protein_index,
				           'databasename'=> $database_basename,
						   'num_pro_per_job'=>$num_protein_per_job,
						   'prot_job_num'=>$prot_job_num,
						   'inclusion_decoy'=>$inclusion_decoy,
						   };
				
		Create_Sort_BashFile($parameters,$dta_path)							
	}
	print "  You submit $num_mass_region jobs for sorting index files\n";
	#Waint For the jobs until they get finished
	Check_Job_stat("sort_",$num_mass_region);
	
	#Merge all the files		
	print "\n  Mergering files\n";
	my $merge_job_num = $num_mass_region-1;
	
	my $cmd = "for i in {0..$merge_job_num} \n do\n cat $mass_index.".'$i'." >> $mass_index\n done\n";
	my $FileName = "$dta_path/merge_mass_index.sh";
	LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_mass_index");		
	$cmd = "for i in {0..$merge_job_num} \n do\n cat $peptide_index.".'$i'." >> $peptide_index\n done\n";
	$FileName = "$dta_path/merge_peptide_index.sh";
	LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_peptide_index");
	
### create a script to sum the numbers #############	
	my @summaryfiles=();
	for($i=0;$i<$merge_job_num;$i++)
	{
		my $sumfile = "$mass_index.${i}.summary";
		push(@summaryfiles,$sumfile);
	}
## generate a database summary file#####	
	my $outfile = $mass_index;
	$outfile =~ s/mdx/sdx/;
	summary_db(\@summaryfiles,$outfile);

####################################################	
	print "  You submit 2 jobs for merging index files\n";
	Check_Job_stat("merge_","2");		

	
	
	#Clean Temporary files
	print "\n  Removing temporary files\n";
	system(qq(rm -rf $database_path/temp_*));
	system(qq(rm -rf $mass_index.*));
	system(qq(rm -rf $peptide_index.*));
	print "  Database creation completed\n";
}		

sub summary_db
{
	my ($files,$outfile)=@_;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
	my $fasta_file = $params->{'database_name'};
	my $enzyme_info = $params->{'enzyme_info'};
	my $digestion = $params->{'digestion'};	
	my $mis_cleavage = $params->{'max_mis_cleavage'};
	my $min_mass = $params->{'min_peptide_mass'};	
	my $max_mass = $params->{'max_peptide_mass'};	
	
	my %hash;
	my $total_unique_peptide = 0;
	foreach my $file(@$files)
	{
		open(FILE,$file);
		while(<FILE>)
		{
			chomp $_;
			my @data = split(/\t/,$_);
			if(!defined($hash{$data[0]}{$data[1]}{$data[2]}))
			{
				$hash{$data[0]}{$data[1]}{$data[2]} = $data[3];
			}
			else
			{
				$hash{$data[0]}{$data[1]}{$data[2]} += $data[3];
			}
			$total_unique_peptide += $data[3];	
		}
	}
	close(FILE);
	open(OUTPUT,">$outfile");
	print OUTPUT "JUMP Database index file\n";
	printf OUTPUT "%4d-%02d-%02d %02d:%02d:%02d\n\n",$year+1900,$mon,$mday,$hour,$min,$sec;	
	print OUTPUT "Database_name = ", $fasta_file,"\n";
	print OUTPUT "enzyme_info = ", $enzyme_info,"\n";
	print OUTPUT "digestion = ", $digestion,"\n";
	print OUTPUT "max_mis_cleavage = ", $mis_cleavage,"\n";
	print OUTPUT "max_modification = ", $params->{'max_modif_num'},"\n";	
	print OUTPUT "min_peptide_mass = ", $min_mass,"\n";
	print OUTPUT "max_peptide_mass = ", $max_mass,"\n\n";
	print OUTPUT "The distribution of unique peptides for each group\n";
	print OUTPUT "miscleavage_num\tmodification_num\tpeptide_num\t\fraction\n";
#	my %cum_mis;
#	my %cum_mod;
	foreach my $tryp (sort {$a<=>$b} keys %hash)
	{
		foreach my $mis_num (sort {$a<=>$b} keys %{$hash{$tryp}})
		{
			foreach my $mod_num (sort {$a<=>$b} keys %{$hash{$tryp}{$mis_num}})
			{
				$cum_mis += $hash{$tryp}{$mis_num}{$mod_num};
				
				print OUTPUT $mis_num,"\t",$mod_num,"\t",$hash{$mis_num}{$mod_num},"\t",$hash{$mis_num}{$mod_num}/$total_unique_peptide,"\n";
			}
		}
	}	
	
=head	
	print OUTPUT "miscleavage_num\tfraction\n";
	my $cum_total_mis=0;
	foreach my $mis_num (sort{$a<=>$b} keys %cum_mis)
	{
		
		$cum_total_mis += $cum_mis{$mis_num}; 
		print OUTPUT "mis",$mis_num," = ",$cum_total_mis / $cum_mis{'0'},"\n";		
	}
	print OUTPUT "modification_num\tfraction\n";
	my $cum_total_mod=0;	
	foreach my $mod_num (sort{$a<=>$b} keys %cum_mod)
	{
		$cum_total_mod += $cum_mod{$mod_num};
		print OUTPUT "mod",$mod_num," = ",$cum_total_mod / $cum_mod{'0'},"\n";		
	}
=cut	
	close(OUTPUT);	
}
###################### window part ##########################
print "  Converting .raw into .mzXML file\n";

my $mzXML = $proc_raw->raw2mzXML();
##################### Linux part ############################
print "  Extracting peaks from .mzXML\n";
my $proc_xml = new Spiders::ProcessingMzXML();
$proc_xml ->set_dta_path($dta_path);

$proc_xml ->set_mzXML_file($mzXML);

################### preprocessing #########################
my (%ms_hash,%msms_hash,@mz_array);

$proc_xml ->set_parameter($params);
$proc_xml->generate_hash_dta(\%ms_hash, \%msms_hash, \@mz_array, $params);



printf "\n  There are %d MS and %d MS/MS in the entire run\n", scalar(keys %{$ms_hash{'surveyhash'}}) , scalar(keys %msms_hash)-scalar(keys %{$ms_hash{'surveyhash'}});

if($params->{'bypass_precursor_preprocessing'} == 0)
{
	print "\n  Decharging for each scans\n";
	my $pip = new Spiders::PIP;
	$pip->set_parameter($params);
	$pip->set_origmz_array(\@mz_array);
	$pip->set_origmsms_hash(\%msms_hash);
	$pip->set_isolation_window($params->{'isolation_window'});
	$pip->set_dta_path($dta_path);	
	my $PIPref = $pip->Calculate_PIP();

	my ($charge_dist,$ppi_dist) = $pip->changeMH_folder($PIPref);
}
elsif($params->{'bypass_precursor_preprocessing'} == 1)
{
	print "  Bypassing decharging for identified MS1 precursor ion\n";
	
	my $charge_assignment = $params->{'manual_charge_assignment'};
	$charge_assignment =~ s/\s*//g;
		
	my @charge_array = split(/\,/,$charge_assignment);

	
	my $pip = new Spiders::PIP;
	$pip->set_parameter($params);
	$pip->set_origmz_array(\@mz_array);
	$pip->set_origmsms_hash(\%msms_hash);
	$pip->set_isolation_window($params->{'isolation_window'});
	$pip->set_dta_path($dta_path);	
	my $PIPref = $pip->Calculate_PIP();
	my @dta = glob("$dta_path/*.dta");		
	my $newdir = "\.${dta_path}.1";
	system(qq(mkdir $newdir));
	

	my $dta_number = scalar @dta;
	foreach my $dtafile (@dta)
	{	
		$dtafile =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;
		open (IN, "<$dtafile");
		my @inarray = <IN>;
		close IN;			
		
		my $specscan = $3;
		my $orig_charge = $5;			
		my $prec_mz =  $PIPref->{$specscan}->{1}->{'mono_mz'};
		foreach my $charge (@charge_array)
		{	
			$dtafile = $2 . "." . $specscan . "." . "1" . "." . $charge . ".dta";
			open (OUT, ">$newdir/$dtafile");
			shift @inarray;
			my $MH = sprintf("%.5f",(($prec_mz-$H)*$charge+$H));
			print OUT "$MH $charge\n";
			print OUT @inarray;
			close(OUT);							
		}		
	}
	system(qq(rm -R $dta_path/ >/dev/null 2>&1));
	system(qq(rm -rf $dta_path >/dev/null 2>&1));	
	system(qq(mv $newdir $dta_path >/dev/null 2>&1));		
}


=head
my $decharge = new Spiders::Decharge();
$decharge->set_parameter($params);
$decharge->set_dta_path($dta_path);	
my @dta = glob("$curr_dir/$dta_path/*.dta");
my %chargehash;
my %mzhash;
foreach my $dta_file (@dta)
{
	my $dta_hash = $decharge->create_dtahash($dta_file,\%msms_hash);
	my $charge= $decharge->decharge($dta_file, $dta_hash, \%msms_hash, \%ms_hash, \@mz_array, \%realcharge);
	
	$dta_file =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;		
	$chargehash{$3}=$charge;
	$mzhash{$3}=$decharge->{'_deisotope_mz_hash'};
}
$decharge->changeMH_folder(\%msms_hash,\%chargehash,\%mzhash);
=cut

###### Build dynamic mass tolerance based on intensity ##################
if($params->{'vary_tolerance'}==1)
{
	my $massaccu = new Spiders::MassAccuracy();
	$massaccu->setmsmshash(\%msms_hash);
	$massaccu->setmasstolerance($params->{'peptide_tolerance'});
	$massaccu->setmasstolerance_frac($params->{'peptide_tolerance_frac'});
	my $dynamic_mass_tolerance = $massaccu->calculate_dynamic_tolerance();
	store($dynamic_mass_tolerance,"$curr_dir/$dta_path/.dynamic_mass_tolerance");
}



########################## Start Searching #######################################
print "  Starting database searching\n";

system(qq(mv sequest.params $dta_path >/dev/null 2>&1));

my $job = new Spiders::Job;
$job->set_library_path($library);
$job->set_dta_path("$curr_dir/$dta_path");
$job->set_pip($PIPref);
$job->create_script();

my @file_array = glob("$curr_dir/$dta_path/*.dta");
#my $dta_num_per_file = $params->{'scans_per_job'};

my $job_num = 200;
my $dta_num_per_file = int($#file_array / $job_num) + 1;
#my $dta_num_per_file = 200/($num_dynamic_mod*2);
#my $job_num=int($#file_array/$dta_num_per_file)+1;
$job_num = $#file_array+1 if($#file_array<$job_num);

for($i=0;$i<$job_num;$i++)
{
	open(JOB,">$curr_dir/$dta_path/job_$i.sh") || die "can not open the job files\n";
	my $dta_file_temp="";
	for($j=0;$j<$dta_num_per_file;$j++)
	{
		if(($i*$dta_num_per_file+$j)<=$#file_array)
		{
			$dta_file_temp .= " $file_array[$i*$dta_num_per_file+$j]";
		}
	}
	if($params->{'Job_Management_System'} eq 'LSF')
	{
		print JOB "#BSUB -P prot\n";
		print JOB "#BSUB -q normal\n";
#		print JOB "#BSUB -M 2000\n";
#		print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
		print JOB "#BSUB -eo $curr_dir/$dta_path/$i.e\n";
		print JOB "#BSUB -oo $curr_dir/$dta_path/$i.o\n";
		print JOB "perl $curr_dir/$dta_path/runsearch_shell.pl -param $curr_dir/$parameter -dta_path $curr_dir/$dta_path $dta_file_temp\n";		
	}
	elsif($params->{'Job_Management_System'} eq 'SGE')
	{
		print JOB "#!/bin/bash\n";
        print JOB "#\$ -N job_$i\n";
        print JOB "#\$ -e $curr_dir/$dta_path/$i.e\n";
        print JOB "#\$ -o $curr_dir/$dta_path/$i.o\n";
		print JOB "perl $curr_dir/$dta_path/runsearch_shell.pl -param $curr_dir/$parameter -dta_path $curr_dir/$dta_path $dta_file_temp\n";	
	}
	close(JOB);
}

######### running jobs ######################### 
my $job_list;
if($params->{'cluster'} eq '1')
{
	if($params->{'Job_Management_System'} eq 'LSF')
	{
		for(my $i=0;$i<$job_num;$i++)
		{
			$command_line = qq(cd $dta_path && bsub <job_$i.sh);
			my $job=qx[$command_line];
			chomp $job;
			my $job_id=0;
			if($job=~/Job \<(\d*)\> is/)
			{
				$job_id=$1;
			}
			$job_list->{$job_id}=1;
		}
	}
	elsif($params->{'Job_Management_System'} eq 'SGE')
	{

		for(my $i=0;$i<$job_num;$i++)
		{
			$command_line = qq(cd $dta_path && qsub -cwd job_$i.sh );
			my $job=qx[$command_line];
			chomp $job;
			my $job_id=0;
			if($job=~/Job \<(\d*)\> is/)
			{
				$job_id=$1;
			}
			$job_list->{$job_id}=1;
			my $count = $i+1;
			print "\r  $count jobs were submitted";
			
		}	
	}
}


	system(qq(cp -rf $parameter "$curr_dir/$dta_path/jump.params" >/dev/null 2>&1));	
	print "\n  You submitted $job_num jobs for database search\n";

	Check_Job_stat("job_",$job_num);
	print "\n  Searching finished\n\n";
#if($params->{search_engine} eq 'JUMP')
#{
	print "  Generating pepXML file\n";
	my $spout=new Spiders::SpoutParser();

	my $outXML = $dta_path . ".pepXML"; 
	my (%seqParaHash, %outInforHash);

	$spout->readSpsearchParam(\%seqParaHash,$parameter);
	$spout->parseSpout(\%seqParaHash, \%outInforHash, "$curr_dir/$dta_path");
	$spout->printPepXML(\%seqParaHash, \%outInforHash, "$curr_dir/$dta_path", $outXML, 5);
	print "\n";
#}
=head
print "\n  Generating summary XLSX report file\n";
$spout->set_MS1scanNumber($ms1N);
$spout->set_MS2scanNumber($ms2N);
$spout->set_totalScan($ms1N+$ms2N);
$spout->set_chargeDistribution($charge_dist);
$spout->set_ppiDistribution($ppi_dist);

$spout->printXLSX(\%seqParaHash, \%outInforHash, "$curr_dir/$dta_path",$topMS2_array);

=cut






			
sub Create_Sort_BashFile
{
	
	my ($Params,$dta_path) = @_;	
	my $FileName = "$curr_dir/$dta_path/sort_db_".$Params->{"range"}.".sh";	
	my $cmd = join(" ","perl $curr_dir/$dta_path/Create_Partial_Idx.pl -m",$Params->{"range"},"-job_num",$Params->{"JobNum"},"-dta_path",$Params->{"dta_path"},"-database_path",$Params->{"database_path"},"-mass_index",$Params->{"mass_index"},"-peptide_index", $Params->{"peptide_index"},"-protein_index",$Params->{"protein_index"},"-databasename",$Params->{"databasename"},"-num_pro_per_job",$Params->{"num_pro_per_job"},"-prot_job_num",$Params->{"prot_job_num"},"-inclusion_decoy",$Params->{"inclusion_decoy"},"\n");
	#my $cmd = "perl Create_Partial_Idx.pl -m ".$Params->{"range"}." -job_num ".$Params->{"JobNum"}." -dta_path ".$Params->{"dta_path"}." -curr_dir ".$Params->{"curr_dir"}." -mass_index ".$Params->{"mass_index"}." -peptide_index ".$Params->{"peptide_index"}." -protein_index ".$Params->{"protein_index"}."\n";
	
	LuchParallelJob($FileName,$cmd,$Params->{"GridType"},"sort_db_".$Params->{"range"});

}


sub LuchParallelJob{
	
	my($FileName,$cmd,$GridType,$outputName)= @_;
	my $mem_free = 8 * $partial_job . "G";
	my $max_mem = 16 * $partial_job . "G";
	open(JOB,">$FileName") || die "can not open $FileName\n";
	if($GridType eq 'LSF')
	{	
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q normal\n";
			print JOB "#BSUB -M 20000\n";
			print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
			
			print JOB "#BSUB -eo $curr_dir/$dta_path/$outputName.e\n";
			print JOB "#BSUB -oo $curr_dir/$dta_path/$outputName.o\n";
			print JOB $cmd;		
			close(JOB);
			system(qq(bsub <$FileName >/dev/null 2>&1));	
	}
	if($GridType eq 'SGE')
	{

		print JOB "#!/bin/bash\n";
		print JOB "#\$ \-S /bin/bash\n";  #In our cluster this line is esential for executing some bash commands such as for
		print JOB "#\$ \-N $outputName\n";
		print JOB "#\$ \-e $curr_dir/$dta_path/$outputName.e\n";
		print JOB "#\$ \-o $curr_dir/$dta_path/$outputName.o\n";
		print JOB $cmd;
		close(JOB);
		system(qq(qsub -cwd -pe mpi 4 -l mem_free=$mem_free,h_vmem=$max_mem $FileName >/dev/null 2>&1));
	}
	if($GridType eq 'PBS')
	{

		print JOB "#!/bin/bash\n";
		print JOB "#PBS -N $outputName\n";
		print JOB "#PBS -e $curr_dir/$dta_path/$outputName.e\n"; 
		print JOB "#PBS -o $curr_dir/$dta_path/$outputName.o"; 			
		print JOB $cmd;
		close(JOB);
		system(qq(qsub -cwd $FileName >/dev/null 2>&1));
	}
	close(JOB);
}

sub MergeFile
{
	my $nbRanges = shift;
	my $cmd = "for i in {0..$nbRanges} do; cat "
}
	
sub Check_Job_stat
{
# test whether the job is finished
	my ($jobs_prefix,$job_num) = @_;
	$job_info=1;
    my ($username) = getpwuid($<);
	my $command_line="";
	my $dot = ".";
	while($job_info)
	{
		sleep(10);

		if($params->{'Job_Management_System'} eq 'LSF')
		{
			$command_line =  "bjobs -u $username";
		}
		elsif($params->{'Job_Management_System'} eq 'SGE')
		{
			$command_line =  "qstat -u $username";
		}
		my $job_status=qx[$command_line];
		#print "$job_status\n";
## end of change
		#my $job_status=qx{$command_line 2>&1};
		
		my @job_status_array=split(/\n/,$job_status);
	
		#Consider only the one that we submitted
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			$command_line =  "bjobs -u $username";	

			my $job_status=qx[$command_line];
			my @job_status_array=split(/\n/,$job_status);
			my $job_number = $job_num - scalar (@job_status_array);
			if(scalar (@job_status_array) == 0)
			{
				print "\r  $job_num jobs finished          ";
			}
			else
			{
				print "\r  $job_number jobs finished          ";
				sleep(5);
			}
			if(scalar(@job_status_array)>0)
			{
				$job_info=1;				
			}
			else
			{
				$job_info=0;		
			}			
		}
		elsif($params->{'Job_Management_System'} eq 'SGE')
		{
		
			@job_status_array = grep(/$jobs_prefix/,@job_status_array); 
			if($jobs_prefix eq "job_")
			{
		#		my $check_command = "ls -f $dta_path\/\*.spout \| wc -l";
		#		my @outfile = glob("$dta_path\/\*.spout");

		#		my $outfile_num=scalar @outfile;
				#chomp $outfile_num;
=head				
				$command_line =  "qstat -u $username | wc -l";	
				my $job_status=qx[$command_line];
				my @job_status_array=split(/\n/,$job_status);
				my $job_number = $job_num - scalar (@job_status_array) + 2;
				
				my $count = $#file_array+1;
				my $outfile_num = $count - ($job_number - 2) * $dta_num_per_file;
				print "\r  Searching $outfile_num out of $count dta files          ";
=cut
				$dot .= ".";
				if(length($dot) == 25)
				{
					$dot = ".";
				}
				print "\r  Searching $dot                                                           ";
			}
			elsif($jobs_prefix eq "job_db_" || $jobs_prefix eq "sort_" || $jobs_prefix eq "merge_")
			{
				if($params->{'Job_Management_System'} eq 'LSF')
				{	
					$command_line =  "bjobs -u $username | grep $jobs_prefix";	
				}
				elsif($params->{'Job_Management_System'} eq 'SGE')
				{
					$command_line =  "qstat -u $username | grep $jobs_prefix";			
				}			
				my $job_status=qx[$command_line];
				my @job_status_array=split(/\n/,$job_status);
				my $job_number = $job_num - scalar (@job_status_array);
				if(scalar (@job_status_array) == 0)
				{
					print "\r  $job_num jobs finished          ";
				}
				else
				{
					print "\r  $job_number jobs finished          ";
				}
				if(scalar(@job_status_array)>0)
				{
					$job_info=1;				
				}
				else
				{
					$job_info=0;		
				}
			}		
		}
		if($job_status=~/No unfinished job found/)
		{
			$job_info=0;
			print "  \n";
		}
		elsif((scalar(@job_status_array))==0)
		{
			$job_info=0;
		}
		elsif($job_status_array[1]=~/PEND/)
		{
			print "\r cluster is busy, please be patient!          ";
			sleep(5);
		}
	}
}


sub check_input
{
	my ($raw_file,$sim_path,$parameter)=@_;

	my $err = new Spiders::Error;
	if(!defined($raw_file))
	{
		$err->no_rawfile_error();
	}
	if(!defined($parameter))
	{
		$parameter = "spSearch.params";
	}
}

sub usage {

print <<"EOF";
	
################################################################
#                                                              #
#       **************************************************     # 
#       ****                                          ****     # 
#       ****                 JUMP                     ****     # 
#       ****              version 1.16                ****     # 
#       ****        Xusheng Wang / Junmin Peng        ****     # 
#       ****         Copyright (C) 2012 - 2013        ****     # 
#       ****            All rights reserved           ****     # 
#       ****                                          ****     # 
#       **************************************************     # 
#                                                              #
################################################################

Usage: $progname <options> rawfile.raw
   -param              input the parameter file
   -rawfile            input the raw file
   -help               print this measly help
EOF
exit 1;
}
