#!/usr/bin/perl  -I /home/xwang4/scripts -I /usr/local/lib/perl5

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
use Spiders::pepXML;
use Spiders::PIP;


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

###################### window part ##########################
print "  Converting .raw into .mzXML file\n";
my $proc_raw = new Spiders::ProcessingRAW();
$proc_raw->set_raw_file($raw_file);
my $dta_path = $proc_raw->get_rawfile_basename();

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



print "\n  Decharging for each scans\n";
my $pip = new Spiders::PIP;
$pip->set_parameter($params);
$pip->set_origmz_array(\@mz_array);
$pip->set_origmsms_hash(\%msms_hash);
$pip->set_isolation_window($params->{'isolation_window'});
$pip->set_dta_path($dta_path);	
my $PIPref = $pip->Calculate_PIP();

$pip->changeMH_folder($PIPref);



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
########################################################################

#my $fasta_file = $params->{'database_name'};
my $databasename = $params->{'database_name'};
my $database_path = dirname($databasename); 
my $database_basename = basename($databasename);

if(!(-e ($databasename) and $databasename=~/.mdx/))
{
	print "  Creating database\n";
	my $mass_index = $databasename . ".mdx";

	my $peptide_index = $databasename . ".pdx";
	my $protein_index = $databasename . ".prdx";
	if(-e($mass_index))
	{
		print "Do you want to remove the old database with same name?\n";
		chomp(my $choice = <STDIN>);
		if($choice eq "yes" || $choice eq "y")
		{
			system(qq(rm -f $mass_index >/dev/null 2>&1));
			system(qq(rm -f $peptide_index >/dev/null 2>&1));
			system(qq(rm -f $protein_index >/dev/null 2>&1));
		}	
	}
	$params->{'database_name'} = $mass_index;
	################# create database using cluster system ########################
	my $num_dynamic_mod = 1;
	foreach my $key (keys %$params)
	{
		if($key=~/dynamic_/)
		{
			$num_dynamic_mod++;
		}
	}
	
	my $num_protein_per_job = int(200000/($num_dynamic_mod));
	my $protein_num=0;
	my $k=0;
	

	open(FASTA, $databasename) || die "can not open the database\n";
	while(<FASTA>)
	{
		$protein_num++ if($_=~/^\>/);
		if(($protein_num % $num_protein_per_job)==1)
		{
			if($_=~/^\>/)
			{
				$k++;
				open(FASTATEMP,">$database_path/temp_${k}_${database_basename}");
			}
			print FASTATEMP "$_";

		}
		else
		{
			print FASTATEMP "$_";		
		}
	}
	close(FASTA);

########### make job file ##############	
	my $job = new Spiders::Job();
	my $abs_parameter = abs_path ($parameter);
	my $num_mass_region = 20 * $num_dynamic_mod;
#	my $num_mass_region = 20;
	$job->make_createdb_script("$curr_dir/$dta_path",$abs_parameter,$num_mass_region);

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
			system(qq(qsub -cwd -pe mpi 16 -l mem_free=4G,h_vmem=6G "$dta_path/job_db_$i.sh" >/dev/null 2>&1));
		}		
	}
	print "  You submit $job_num jobs for creating index files\n";
	
	Check_Job_stat("job_db_",$job_num);

###### Merge all temp_ files into big data file #######################
	
	my $j=0;
	my $l=0;
	my $prev_run = 0;
	print "\n  Sorting indexes\n";

	$job->make_partialidx_script("$curr_dir/$dta_path");

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
						   'num_pro_per_job'=>$num_protein_per_job
						   } ;
				
		Create_Sort_BashFile($parameters)							
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
	print "  You submit 2 jobs for merging index files\n";
	Check_Job_stat("merge_","2");		
	
	#Clean Temporary files
	system(qq(rm -rf $database_path/temp_*));
	system(qq(rm -rf $mass_index.*));
	system(qq(rm -rf $peptide_index.*));	
}		


########################## Start Searching #######################################
print "  Starting database searching\n";

system(qq(cp sequest.params $dta_path));
my $job = new Spiders::Job;
$job->set_dta_path("$curr_dir/$dta_path");
$job->set_pip($PIPref);
$job->create_script();

my @file_array = glob("$curr_dir/$dta_path/*.dta");
my $job_num=int($#file_array/100)+1;

for($i=0;$i<$job_num;$i++)
{
	open(JOB,">$curr_dir/$dta_path/job_$i.sh") || die "can not open the job files\n";
	my $dta_file_temp="";
	for($j=0;$j<100;$j++)
	{
		if(($i*100+$j)<=$#file_array)
		{
			$dta_file_temp .= " $file_array[$i*100+$j]";
		}
	}
	if($params->{'Job_Management_System'} eq 'LSF')
	{
		print JOB "#BSUB -P prot\n";
		print JOB "#BSUB -q normal\n";
		print JOB "#BSUB -M 20000\n";
		print JOB "#BSUB -R \"rusage[mem=20000]\"\n";			
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
			$job_list->{$job}=1;
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
			$job_list->{$job}=1;
		}	
	}
}

Check_Job_stat("job_",$job_num);
print "  Searching finished\n\n";

=head
my $ss = new Spiders::SummaryResults();

my $mainhash=$ss->parse_spOut_files_new("$curr_dir/$dta_path/");
my $outhash=$ss->parse_outfile("$curr_dir/$dta_path/");

for my $spout (keys %$mainhash)
{
		print $spout,"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"MH"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"peptide_E_value"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"reference"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"peptide"},"\t"; 
		print $out,"\t",$outhash->{$spout}->{'XCorr'},"\t",$outhash->{$spout}->{'dCn'},"\t",$outhash->{$spout}->{'peptide'},"\n"; 
		
}
open(RESULT1,">XcorrEscore.txt");
$ss->compare_psearch_sequest_scores();


system(qq(R --save <fig1.r));
system(qq(R --save <fig2.r));
system(qq(R --save <fig3.r));
=cut


################
=head
print "  Generating pepXML file\n";
my $pepxml = new Spiders::pepXML();

my $outXML = $dta_path . ".pepXML"; 
my (%seqParaHash, %outInforHash);

$pepxml->readSpsearchParam(\%seqParaHash,$parameter);
$pepxml->parseSpout(\%seqParaHash, \%outInforHash, "$curr_dir/$dta_path");
$pepxml->printPepXML(\%seqParaHash, \%outInforHash, "$curr_dir/$dta_path", 1, $outXML);

############## Summarize the data #####################
my $ss = new Spiders::SummaryResults('-path'=>"$curr_dir/$dta_path/");

open(SUM,">Summary.csv");
my $mainhash=$ss->parse_spOut_files_v4("$curr_dir/$dta_path/");
#my $outhash=$ss->parse_outfile("Rat_B_100ng_Q");
foreach (keys %$mainhash)
{

        print SUM $_,"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"reference"},$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"peptide_E_value"},"\t",$mainhash->{$_}->{"outhash"}->{"Order"}->{"1"}->{"Weight_E_value"},"\n";
}


my $comp_score = $ss->compare_psearch_sequest_scores();
open(RESULT1,">XcorrEscore.txt");
print RESULT1 "Xcorr\tPeptide_E_score\tMatch\n";
foreach my $file (keys %$comp_score)
{
        if(!defined($comp_score->{$file}->{"Escore"}))
        {
                $comp_score->{$file}->{"Escore"}=0;
        }
        if(!defined($comp_score->{$file}->{"Xcorr"}))
        {
                $comp_score->{$file}->{"Xcorr"}=0;
        }
        print RESULT1 $file,"\t",$comp_score->{$file}->{"Xcorr"},"\t",$comp_score->{$file}->{"Escore"},"\t",$comp_score->{$file}->{"match"},"\n";
}

#system(qq(R --save <fig3.r));
close(RESULT1);

open(RESULT2,">WeightedDist.txt");
print RESULT2 "Evalue_Target\tNumber_of_IDs\tXcorr_Decoy\tNumber_of_IDs";
my $score = $ss->score_distribution($mainhash);

my %merge_score;
foreach (keys %{$score->{"Decoy"}})
{
        $merge_score{$_}=1;
}

foreach (keys %{$score->{"Target"}})
{
        $merge_score{$_}=1;
}

open(FDRRESULT,">FDR_Escore.txt");
print FDRRESULT "E_score\tFDR\n";
my $sum_target = 0;
my $sum_decoy = 0;
foreach (reverse sort {$a<=>$b} keys %merge_score)
{
        $score->{"Target"}->{$_}=0 if(!defined($score->{"Target"}->{$_}));
        $score->{"Decoy"}->{$_}=0 if(!defined($score->{"Decoy"}->{$_}));
        next if ($_==0);
        if($_<=100)
        {
                my $FDR = $score->{"Decoy"}->{$_}/$score->{"Target"}->{$_};
                print FDRRESULT $_,"\t",$FDR,"\n";
        }
        $sum_target = $score->{"Target"}->{$_}+$sum_target;
        $sum_decoy = $score->{"Decoy"}->{$_}+$sum_decoy;
        print RESULT2 $_,"\t",$sum_target,"\t",$sum_decoy,"\n";
}
=cut
########################


				
sub Create_Sort_BashFile
{
	
	my $Params = shift;	
	my $FileName = "$curr_dir/$dta_path/sort_db_".$Params->{"range"}.".sh";	
	my $cmd = join(" ","perl $curr_dir/$dta_path/Create_Partial_Idx.pl -m",$Params->{"range"},"-job_num",$Params->{"JobNum"},"-dta_path",$Params->{"dta_path"},"-database_path",$Params->{"database_path"},"-mass_index",$Params->{"mass_index"},"-peptide_index", $Params->{"peptide_index"},"-protein_index",$Params->{"protein_index"},"-databasename",$Params->{"databasename"},"-num_pro_per_job",$Params->{"num_pro_per_job"},"\n");
	#my $cmd = "perl Create_Partial_Idx.pl -m ".$Params->{"range"}." -job_num ".$Params->{"JobNum"}." -dta_path ".$Params->{"dta_path"}." -curr_dir ".$Params->{"curr_dir"}." -mass_index ".$Params->{"mass_index"}." -peptide_index ".$Params->{"peptide_index"}." -protein_index ".$Params->{"protein_index"}."\n";
	
	LuchParallelJob($FileName,$cmd,$Params->{"GridType"},"sort_db_".$Params->{"range"});

}


sub LuchParallelJob{
	
	my($FileName,$cmd,$GridType,$outputName)= @_;
	
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
		system(qq(qsub -cwd -cwd -pe mpi 16 -l mem_free=4G,h_vmem=6G $FileName >/dev/null 2>&1));
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
			my $job_number = $job_num - scalar (@job_status_array) + 1;
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
				my $check_command = "ls -f $dta_path\/\*.spout \| wc -l";
				my @outfile = glob("$dta_path\/\*.spout");

				my $outfile_num=scalar @outfile;
				#chomp $outfile_num;
				print "\r  Searching $outfile_num out of $#file_array dta files          ";
			}
			elsif($jobs_prefix eq "job_db_" || $jobs_prefix eq "sort_" || $jobs_prefix eq "merge_")
			{
				if($params->{'Job_Management_System'} eq 'LSF')
				{	
					$command_line =  "bjobs -u $username";	
				}
				elsif($params->{'Job_Management_System'} eq 'SGE')
				{
					$command_line =  "qstat -u $username";			
				}			
				my $job_status=qx[$command_line];
				my @job_status_array=split(/\n/,$job_status);
				my $job_number = $job_num - scalar (@job_status_array) + 2;
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
