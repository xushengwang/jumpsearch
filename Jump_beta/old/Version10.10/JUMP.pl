#!/usr/bin/perl  -I /home/xwang4/scripts -I /usr/local/lib/perl5

################################################################
#                                                              #
#       **************************************************     # 
#       **** Multiplexing program                     ****     # 
#       ****                                          ****     # 
#       ****Copyright (C) 20212 - Xusheng Wang        ****     # 
#       ****all rights reserved.                      ****     # 
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
use Spiders::SummaryResults;


my $VERSION = 1.01;

my $progname = $0;
# remove path from our name
$progname =~ s@(.*)/@@i;



my ($help,$parameter);
GetOptions('-help|h'=>\$help,
			'-param=s'=>\$parameter,
		);


######## predefine the input #####################
my $raw_file=$ARGV[0];
usage() if ($help);
check_input($raw_file,\$parameter);
######### programming starting information #############
print "\nInitializing JUMP program\n";

###### Get working directory #####################
my $curr_dir = getcwd;

my $p = Spiders::Params->new('-path'=>$parameter);
my $params=$p->parse_param();
###################### window part ##########################
print "Converting .raw into .mzXML file\n";
my $proc_raw = new Spiders::ProcessingRAW();
$proc_raw->set_raw_file($raw_file);
my $dta_path = $proc_raw->get_rawfile_basename();

my $mzXML = $proc_raw->raw2mzXML();
##################### Linux part ############################
print "Extracting peaks from .mzXML\n";
my $proc_xml = new Spiders::ProcessingMzXML();
$proc_xml ->set_dta_path($dta_path);

$proc_xml ->set_mzXML_file($mzXML);

################### preprocessing #########################
print "Gathering scan hashes ... \n";

my (%ms_hash,%msms_hash,@mz_array);

$proc_xml ->set_parameter($params);
$proc_xml->generate_hash_dta(\%ms_hash, \%msms_hash, \@mz_array, $params);

printf "\nThere are %d MS and %d MS/MS in the entire run\n", scalar(keys %{$ms_hash{'surveyhash'}}), scalar(keys %msms_hash);

print "  Decharging for each scans\n";
my $decharge = new Spiders::Decharge();
$decharge->set_parameter($params);
$decharge->set_dta_path($dta_path);	
my @dta = glob("$curr_dir/$dta_path/*.dta");
my %chargehash;
foreach my $dta_file (@dta)
{
	my $dta_hash = $decharge->create_dtahash($dta_file,\%msms_hash);
	my $charge= $decharge->decharge($dta_file, $dta_hash, \%msms_hash, \%ms_hash, \@mz_array, \%realcharge);
	$dta_file =~ s/(([A-Za-z0-9\_\-]+)\.(\d+)\.(\d+)\.(\d+).dta)\Z/$1/;		
	$chargehash{$3}=$charge;
}
$decharge->changeMH_folder(\%msms_hash,\%chargehash);


my $fasta_file = $params->{'database_name'};
my $databasename = basename($fasta_file);
if(!-e ("$curr_dir/$dta_path/.database/$databasename.mdx"))
{

	print "Creating database\n";
	
	my $mass_index = $dta_path . "/.database/$databasename.mdx";
	my $peptide_index = $dta_path . "/.database/$databasename.pdx";
	my $protein_index = $dta_path . "/.database/$databasename.prdx";
	
	if (!(-e "$dta_path/.database"))
	{
		system(qq(mkdir $curr_dir/$dta_path/.database));
	}
	else
	{
		system(qq(rm -rf $curr_dir/$dta_path/.database/*));
	}

	################# create database using cluster system ########################
	my $num_protein_per_job = 1000;
	my $protein_num=0;
	my $k=0;
	open(FASTA, $fasta_file) || die "can not open the database\n";
	while(<FASTA>)
	{
		$protein_num++ if($_=~/^\>/);
		if(($protein_num % $num_protein_per_job)==1)
		{
			if($_=~/^\>/)
			{
				$k++;
				open(FASTATEMP,">$curr_dir/$dta_path/.database/temp_${k}_$databasename");
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
	$job->make_createdb_script("$curr_dir/$dta_path",$abs_parameter);

######### submit job file ##############	
	my $job_num=int($protein_num/$num_protein_per_job)+1;	
	for($i=1;$i<=$job_num;$i++)
	{
		open(JOB,">$curr_dir/$dta_path/job_db_$i.sh") || die "can not open the job db files\n";
		if($params->{'Job_Management_System'} eq 'LSF')
		{
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q priority\n";
			print JOB "#BSUB -eo $curr_dir/$dta_path/$i.e\n";
			print JOB "#BSUB -oo $curr_dir/$dta_path/$i.o\n";
			print JOB "perl $curr_dir/$dta_path/create_db.pl $curr_dir/$dta_path/.database/temp_${i}_$databasename\n";		
			close(JOB);
			system(qq(cd $dta_path && bsub <job_db_$i.sh));	
		}
		if($params->{'Job_Management_System'} eq 'SGE')
		{
			print JOB "#!/bin/bash\n";
			print JOB "#\$ \-N job_$i\n";
			print JOB "#\$ \-e $curr_dir/$dta_path/$i.e\n";
			print JOB "#\$ \-o $curr_dir/$dta_path/$i.o\n";
			print JOB "perl $curr_dir/$dta_path/create_db.pl $curr_dir/$dta_path/.database/temp_${i}_$databasename\n";
			close(JOB);
			system(qq(qsub -cwd "$dta_path/job_db_$i.sh"));
		}
		
	}
	Check_Job_stat("job_");
	
###### Merge all temp_ files into big data file #######################

	my $j=0;
	my $l=0;
	my $prev_run = 0;
	print "Sorting indexes\n";

	$job->make_partialidx_script("$curr_dir/$dta_path");
	
	for($m=0;$m<20;$m++)
	{ 		
		#print $m,"\n";
		my ($mass_hash,$peptide_hash,$protein_hash);	
		#create bash files
		my $parameters = {'range'=> $m,				
						  'GridType'=> $params->{'Job_Management_System'},
				          'JobNum'=>$job_num,
				          'curr_dir'=> $curr_dir,
				          'dta_path'=> $dta_path,
				          'mass_index'=>  $mass_index.".$m",
				          'peptide_index'=> $peptide_index.".$m",
				          'protein_index'=> $protein_index,
				           'databasename'=>$databasename,
						   'num_pro_per_job'=>$num_protein_per_job
						   } ;
				
		Create_Sort_BashFile($parameters)							
	}
	
	#Waint For the jobs until they get finished
	Check_Job_stat("sort_");
	
	#Merge all the files		
	print "Mergering files\n";
	
	my $cmd = "for i in {0..19} \n do\n cat $mass_index.".'$i'." >> $mass_index\n done\n";
	my $FileName = "$dta_path/merge_mass_index.sh";
	LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_mass_index");		
	
	$cmd = "for i in {0..19} \n do\n cat $peptide_index.".'$i'." >> $peptide_index\n done\n";
	$FileName = "$dta_path/merge_peptide_index.sh";
	LuchParallelJob($FileName,$cmd,$params->{'Job_Management_System'},"merge_peptide_index");
		
	#Clean Temporary files
	#system(qq(rm -rf $curr_dir/$dta_path/.database/temp_*));
	#system(qq(rm -rf $mass_index.*));
	#system(qq(rm -rf $peptide_index.*));
	
}		
Check_Job_stat("merge_");

########################## Start Searching #######################################
print "Starting database searching\n";

system(qq(cp sequest.params $dta_path));
my $job = new Spiders::Job;
$job->set_dta_path("$curr_dir/$dta_path");
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
		print JOB "#BSUB -q priority\n";
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
if($params->{'cluster'} eq '1')
{
	if($params->{'Job_Management_System'} eq 'LSF')
	{
		my $job_list;
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
		my $job_list;
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

Check_Job_stat("job_");



my $ss = new Spiders::SummaryResults();

my $mainhash=$ss->parse_spOut_files_new("$curr_dir/$dta_path/");
my $outhash=$ss->parse_outfile("$curr_dir/$dta_path/");

for my $spout (keys %$mainhash)
{
		print $spout,"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"MH"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"peptide_E_value"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"reference"},"\t",$mainhash->{$spout}->{"outhash"}->{"Order"}->{1}->{"peptide"},"\t"; 
		print $out,"\t",$outhash->{$spout}->{'XCorr'},"\t",$outhash->{$spout}->{'dCn'},"\t",$outhash->{$spout}->{'peptide'},"\n"; 
}

				
sub Create_Sort_BashFile
{
	
	my $Params = shift;	
	my $FileName = "$curr_dir/$dta_path/sort_db_".$Params->{"range"}.".sh";	
	my $cmd = join(" ","perl $curr_dir/$dta_path/Create_Partial_Idx.pl -m",$Params->{"range"},"-job_num",$Params->{"JobNum"},"-dta_path",$Params->{"dta_path"},"-curr_dir",$Params->{"curr_dir"},"-mass_index",$Params->{"mass_index"},"-peptide_index", $Params->{"peptide_index"},"-protein_index",$Params->{"protein_index"},"-databasename",$Params->{"databasename"},"-num_pro_per_job",$Params->{"num_pro_per_job"},"\n");
	#my $cmd = "perl Create_Partial_Idx.pl -m ".$Params->{"range"}." -job_num ".$Params->{"JobNum"}." -dta_path ".$Params->{"dta_path"}." -curr_dir ".$Params->{"curr_dir"}." -mass_index ".$Params->{"mass_index"}." -peptide_index ".$Params->{"peptide_index"}." -protein_index ".$Params->{"protein_index"}."\n";
	
	LuchParallelJob($FileName,$cmd,$Params->{"GridType"},"sort_db_".$Params->{"range"});

}


sub LuchParallelJob{
	
	my($FileName,$cmd,$GridType,$outputName)= @_;
	
	open(JOB,">$FileName") || die "can not open $FileName\n";
	if($GridType eq 'LSF')
	{	
			print JOB "#BSUB -P prot\n";
			print JOB "#BSUB -q priority\n";
			print JOB "#BSUB -eo $curr_dir/$dta_path/$outputName.e\n";
			print JOB "#BSUB -oo $curr_dir/$dta_path/$outputName.o\n";
			print JOB $cmd;		
			close(JOB);
			system(qq(bsub < $FileName));	
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
		system(qq(qsub -cwd $FileName));
	}
	if($GridType eq 'PBS')
	{

		print JOB "#!/bin/bash\n";
		print JOB "#PBS -N $outputName\n";
		print JOB "#PBS -e $curr_dir/$dta_path/$outputName.e\n"; 
		print JOB "#PBS -o $curr_dir/$dta_path/$outputName.o"; 			
		print JOB $cmd;
		close(JOB);
		system(qq(qsub -cwd $FileName));
		
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
	my $jobs_prefix = shift;
	$job_info=1;
    my ($username) = getpwuid($<);

	while($job_info)
	{
		sleep(10);
		#my $command_line =  "bjobs -u $username";
		$command_line =  "qstat -u $username";
		my $job_status=qx[$command_line];
		#print "$job_status\n";
## end of change
		#my $job_status=qx{$command_line 2>&1};
		
		my @job_status_array=split(/\n/,$job_status);
	
		#Consider only the one that we submitted
		@job_status_array = grep(/$jobs_prefix/,@job_status_array); 
		
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
			print "\r cluster is busy, please be patient!";
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
Multiplexed Mass spectrometry program, Version: $VERSION
Usage: $progname <options> rawfile.raw
   -param              input the parameter file
   -help               print this measly help
   -verbose            lots of messages
EOF
exit 1;
}