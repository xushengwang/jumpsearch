#!/usr/bin/perl -w

package Spiders::pepXML;

use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);

$VERSION     = 1.00;
@ISA	 = qw(Exporter);
@EXPORT      = ();

=head Example 
my $pepxml = Spiders::pepXML();

my $outXML = $dta_path . ".pepXML"; 
my (%seqParaHash, %outInforHash);

$pepxml->readSpsearchParam(\%seqParaHash,$parameter);
$pepxml->parseSpout(\%seqParaHash, \%outInforHash, "$curr_dir/$dta_path");
$pepxml->printPepXML(\%seqParaHash, \%outInforHash, "$curr_dir/$dta_path", 1, $outXML);
=cut

sub new{
	my ($class,%arg)=@_;
    my $self = {

    };
    bless $self, $class;
	return $self;
}

sub AA_mass{
	my ($self,$AA)=@_;
	my (%AA_mass);
	$AA_mass{'G'} = 57.021464;    
	$AA_mass{'D'} = 115.02694;
	$AA_mass{'A'} = 71.037114;
	$AA_mass{'Q'} = 128.05858;
	$AA_mass{'S'} = 87.032029;
	$AA_mass{'K'} = 128.09496;
	$AA_mass{'P'} = 97.052764;
	$AA_mass{'E'} = 129.04259;
	$AA_mass{'V'} = 99.068414;
	$AA_mass{'M'} = 131.04048;
	$AA_mass{'T'} = 101.04768;
	$AA_mass{'H'} = 137.05891;
	$AA_mass{'C'} = 103.00919;
	$AA_mass{'F'} = 147.06841;
	$AA_mass{'I'} = 113.08406;                      
	$AA_mass{'L'} = 113.08406;
	$AA_mass{'J'} = 113.08406;
	$AA_mass{'R'} = 156.10111;
	$AA_mass{'N'} = 114.04293;    
	$AA_mass{'Y'} = 163.06333;
	$AA_mass{'W'} = 186.07931;
	return $AA_mass{$AA} ;
}

sub enzymeStat{
	my ($self,$AA, $status)=@_;
	my (%enzymeStat);
	$enzymeStat{'trypsin'}{'cut'}='KR';
	$enzymeStat{'trypsin'}{'nocuts'}='P';
	$enzymeStat{'trypsin'}{'sense'}='C';

	$enzymeStat{'tryptic'}{'cut'}='KR';
	$enzymeStat{'tryptic'}{'nocuts'}='P';
	$enzymeStat{'tryptic'}{'sense'}='C';

	$enzymeStat{'argc'}{'cut'}='R';
	$enzymeStat{'argc'}{'nocuts'}='P';
	$enzymeStat{'argc'}{'sense'}='C';

	$enzymeStat{'aspn'}{'cut'}='D';
	#$enzymeStat{'aspn'}{'nocuts'}='';
	$enzymeStat{'aspn'}{'sense'}='N';

	$enzymeStat{'chymotrypsin'}{'cut'}='YWFM';
	$enzymeStat{'chymotrypsin'}{'nocuts'}='P';
	$enzymeStat{'chymotrypsin'}{'sense'}='C';

	$enzymeStat{'clostripain'}{'cut'}='R';
	$enzymeStat{'clostripain'}{'nocuts'}='-';
	$enzymeStat{'clostripain'}{'sense'}='C';

	$enzymeStat{'cnbr'}{'cut'}='M';
	$enzymeStat{'cnbr'}{'nocuts'}='P';
	$enzymeStat{'cnbr'}{'sense'}='C';

	$enzymeStat{'elastase'}{'cut'}='GVLIA';
	$enzymeStat{'elastase'}{'nocuts'}='P';
	$enzymeStat{'elastase'}{'sense'}='C';

	if (exists($enzymeStat{$AA}{$status})){	return $enzymeStat{$AA}{$status};}
	else { return ''; }
}

sub gettime{
   my ($self,$sec, $min, $hour, $day, $mon, $year) = localtime(time);
   my $now_date=join("-",($year+1900,$mon+1,$day));
   my $now_time=join(":",($hour,$min,$sec));
   return ($now_date, $now_time);
}

sub parseSpout
{
	my ($self,$seqParaHash, $outInforHash, $run)=@_;

	#deal with out files
	opendir(DIR,$run); my @out = grep {/\.spout\Z/} readdir(DIR);
	closedir(DIR);

	if (scalar(@out) == 0){ die "     THERE ARE NO .spout FILES IN DIRECTORY $run !!!!\n";}
	my ($orignum, $count) = (scalar(@out), 0);

	#extract general params from an out file
	$self->genParOut($seqParaHash,$run,$out[0]);


	#read outfiles
	for my $outfile (@out)
	{
		$count++;
		print "\rGathering information from $count of $orignum outfiles     ";
	
		$self->readOutfile($seqParaHash, $outInforHash, $run, $outfile);
	}
}

sub readSpsearchParam
{
	my ($self,$seqPara,$paramF)=@_;
	open(SQT,"$paramF") or die 'Cannot open $paramF\n';
	
	my @array1=('peptide_mass_tolerance','fragment_ion_tolerance','ion_series','max_num_differential_AA_per_mod','num_output_lines','remove_precursor_peak','ion_cutoff_percentage',
	'match_peak_count','match_peak_allowed_error','match_peak_tolerance','protein_mass_filter');
	$seqPara->{'optionalParams'}=\@array1;
	
	my @array2=('xcorr','deltacn','deltacnstar','spscore','sprank','lSideMass',   'rSideMass',  'PeptideEvalue', 'TagSeq',  'WeightedEvalue',  'TagsNum');
	$seqPara->{'hitOptionalParams'}=\@array2;
	
	my ($DB_local_path);
	while (<SQT>){
		#last if (/SEQUEST_ENZYME_INFO/);
		next if (/^#/ || /^\s+\Z/);
		chomp; 	s/[;#].*//; s/\s+\Z//;		#Q: remove empty lines
		my ($key, $value) = split(" = ", $_);
		next if (!defined($value));
		
		#push(@$seqPara{'optionalParams'},$key);
		$$seqPara{$key}=$value;
	}
	
	#Tryptic KR P
	$$seqPara{'enzyme_info'} =~ /^(\w+)\s(\w+)/; $$seqPara{'enzyme_info'} = $1;
	$$seqPara{'enzyme'}=lc($$seqPara{'enzyme_info'});
	$$seqPara{'misClv'}=$$seqPara{'max_mis_cleavage'};
	$$seqPara{'frmTol'}=$$seqPara{'frag_mass_tolerance'};
	
	close(SQT);
}

sub genParOut
{
	my ($self,$paraHash,$run,$outfile)=@_;
	
	my $fullpath = "$run\/$outfile"; 		#can have a problem
	open (IN, $fullpath);  
	
	
	my $mod;
	while (<IN>)
	{
		chomp;
		if (/^Database\s=/)
		{
			#Database = /home/xwang4/JUMP_database/ratmouse/ratmouse_con.fasta.mdx
			my @tmp=split(/\s=\s/,$_);
			$$paraHash{'first_database_name'}=$tmp[1];
		}
		elsif (/^ion series ABCDVWXYZ/)
		{
			my @tmp=split(/\:/,$_);
			$$paraHash{'ion_series'}=$tmp[1];
			
			$mod=<IN>;
			chomp($mod);
		}
	}
		
	#$$paraHash{'enzyme'}='';
	#$$paraHash{'misClv'}='';
	#$$paraHash{'frmTol'}='';
	
	#modification
	if (defined($mod))
	{
		$mod =~ s/\s*\/\s*/\//g; $mod =~ s/(\s)\s+/$1/g; $mod =~ s/^\s+//; $mod =~ s/\s+\Z//;
		#dynamic_M@=15.99492  add_C_Cysteine=57.02146  
		my @tmp=split(/\s/,$mod);
		
		for my $term (@tmp)
		{
			if ($term =~ /dynamic_(\w+)([\@\#\*\^\~\$])\=([0-9\.]+)/)	#dynamic modification
			{
				my ($AA,$symble,$massdiff)=($1,$2,$3);
		
				$$paraHash{'mod'}=1;
			
				my @tmp2=split(//,$AA);
			
				for (my $i=0; $i<scalar(@tmp2); $i++)
				{
					if (defined($$paraHash{'modification'}{'dynamic'}{'count'})) {$$paraHash{'modification'}{'dynamic'}{'count'}++;}
					else {$$paraHash{'modification'}{'dynamic'}{'count'}=1;}
				
					$$paraHash{'modification'}{'dynamic'}{'AA'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=$tmp2[$i];
					$$paraHash{'modification'}{'dynamic'}{'symble'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=$symble;
					$$paraHash{'modification'}{'dynamic'}{'massdiff'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=$massdiff;
					$$paraHash{'modification'}{'dynamic'}{'mass'}[$$paraHash{'modification'}{'dynamic'}{'count'}]=$self->AA_mass($tmp2[$i])+$massdiff;
				}
			}
			elsif ($term =~ /add_(\w)_(\w+)\=([0-9\.]+)/)					#static modification
			{			
				$$paraHash{'mod'}=1;
				if (defined($$paraHash{'modification'}{'static'}{'count'})) {$$paraHash{'modification'}{'static'}{'count'}++;}
				else {$$paraHash{'modification'}{'static'}{'count'}=1;}
			
				$$paraHash{'modification'}{'static'}{'AA'}[$$paraHash{'modification'}{'static'}{'count'}]=$1;
				$$paraHash{'modification'}{'static'}{'massdiff'}[$$paraHash{'modification'}{'static'}{'count'}]=$3;
				$$paraHash{'modification'}{'static'}{'mass'}[$$paraHash{'modification'}{'static'}{'count'}]=$3+$self->AA_mass($1);
			}
		}
		
	}
	else {print "No modification infor provided.\n";}
	
	close IN;
	
}

sub readOutfile
{
	my ($self,$paraHash,$outInfor, $run, $outfile)=@_;
	my $maxHitCondiered=5;			#max number of hits to be considered
	
	my $fullpath = "$run\/$outfile"; 		#can have a problem
	open (IN, $fullpath);  
	#Rat_B_100ng_Q.10000.10000.3.spout
	$outfile =~ /\.(\d+)\.(\d)\.spout/;  
	my ($scan,$charge) = ($1,$2);
	$$outInfor{$outfile}{'scan'}=$scan;
	$$outInfor{$outfile}{'charge'}=$charge;
	#print "$scan,$charge\n";
	while (<IN>)
	{
		if (/^Precursor mass\s*=\s*([0-9\.]+)/)
		{
			$$outInfor{$outfile}{'MHmass'}=$1;
		}
		elsif (/^-+\n/)
		{
			last;	#here comes the hits
		}
	}
	unless (defined($$outInfor{$outfile}{'MHmass'})) {die "MHmass not found in $fullpath\n";}
	
	my $i=0; my $count;	
	$$outInfor{$outfile}{'hitShown'}=0;
	while (<IN>)
	{
		last if (/^\s*$/);					#?
		last if (/All identified peptide/);
		chomp; #Remove return line
		s/^\s+//; #Remove leading spaces
		my @vals = split /\s+/,$_;
		
		if(scalar(@vals)>10){       die "Error parsing file $outfile at line:  Some values are missing\n"; }
		elsif (scalar(@vals)==10)		#hit lines
		{
			$i++; $$outInfor{$outfile}{'hitShown'}++; #hit counts
			$count=0;	#for alternative proteins
			#Order   (M+H)+   lSideMass   rSideMass  PeptideEvalue TagSeq  WeightedEvalue  TagsNum      Reference                       Peptide    
			#-------------------------------------------------------------------------------------------------------------------------------------------------
			#1      2127.994    983.4356    1145.5659    0.6093    S           1.89           1        ##Decoy__sp|P51650-2|SSDH_RAT   K.ELYEDIGYSKGGYCVYVL.-        
			my $peptide = $vals[9];
			#K.ELYEDIGYSKGGYCVYVL.-
            $peptide =~ s/([A-Z\-])\.([A-Z\@\#\*\^\~\$]+)\.([A-Z\-])/$2/; 
			$$outInfor{$outfile}{'prAA'}[$vals[0]]=$1;
			$$outInfor{$outfile}{'nxAA'}[$vals[0]]=$3;
			
			if ( defined($$outInfor{$outfile}{'fullPeptide'}[$i-1]) and $peptide eq $$outInfor{$outfile}{'fullPeptide'}[$i-1]) #same peptide
			{
				$i--; $$outInfor{$outfile}{'hitShown'}--; #hit counts
				next;
			}
			
			#print "\n$outfile\n$peptide,$$outInfor{$outfile}{'prAA'}[$i],$$outInfor{$outfile}{'nxAA'}[$i]\n";
			
			$$outInfor{$outfile}{'protein'}[$vals[0]]=$vals[8];
			
			
			$$outInfor{$outfile}{'hitRank'}[$vals[0]]=$vals[0];
			$$outInfor{$outfile}{'rank'}[$vals[0]]=$vals[0];
			#$$outInfor{$outfile}{'sprank'}[$vals[0]]='';
			$$outInfor{$outfile}{'expmass'}[$vals[0]]=$vals[1];
			#$$outInfor{$outfile}{'xcorr'}[$vals[0]]='';
			
			$$outInfor{$outfile}{'lSideMass'}[$vals[0]]=$vals[2];
			$$outInfor{$outfile}{'rSideMass'}[$vals[0]]=$vals[3];
			$$outInfor{$outfile}{'PeptideEvalue'}[$vals[0]]=$vals[4];
			$$outInfor{$outfile}{'TagSeq'}[$vals[0]]=$vals[5];
			$$outInfor{$outfile}{'WeightedEvalue'}[$vals[0]]=$vals[6];
			$$outInfor{$outfile}{'TagsNum'}[$vals[0]]=$vals[7];
			$$outInfor{$outfile}{'fullPeptide'}[$i]=$peptide;
			
			
			
			$$outInfor{$outfile}{'deltacn'}[$vals[0]]=0;
			#$$outInfor{$outfile}{'deltacnstar'}[$vals[0]]=0;
			if ($vals[0]>1)
			{
				if ($$outInfor{$outfile}{'WeightedEvalue'}[$vals[0]-1] !=0 )
				{
					$$outInfor{$outfile}{'deltacn'}[$vals[0]-1]=($$outInfor{$outfile}{'WeightedEvalue'}[$vals[0]-1]-$$outInfor{$outfile}{'WeightedEvalue'}[$vals[0]])/$$outInfor{$outfile}{'WeightedEvalue'}[$vals[0]-1];
				}
			}
			
			
			$$outInfor{$outfile}{'matchIons'}[$vals[0]]=0;
			$$outInfor{$outfile}{'totalIons'}[$vals[0]]=0;
			$$outInfor{$outfile}{'red'}[$vals[0]]=0;
			
			
			
		#consider modifications:
		if (defined($$paraHash{'modification'}{'dynamic'}{'count'}))		#dynamic mod:
		{
			$$outInfor{$outfile}{'dynamicMod'}{'found'}[$vals[0]]=0;
			my @tmp=split(//,$peptide);
			
			my $j=0;
			while ($j<scalar(@tmp))
			#for (my $j=0; $j<scalar(@tmp); $j++)
			{
				for (my $k=1; $k<=$$paraHash{'modification'}{'dynamic'}{'count'}; $k++)
				{
					if ($tmp[$j] eq $$paraHash{'modification'}{'dynamic'}{'AA'}[$k])
					{
						if ($j+1<scalar(@tmp) and $tmp[$j+1] eq $$paraHash{'modification'}{'dynamic'}{'symble'}[$k])
						{	
							$j++;
							$$outInfor{$outfile}{'dynamicMod'}{'found'}[$vals[0]]++;
							$$outInfor{$outfile}{'dynamicMod'}{'modIndex'}[$vals[0]][$$outInfor{$outfile}{'dynamicMod'}{'found'}[$vals[0]]]=$k;
							$$outInfor{$outfile}{'dynamicMod'}{'pos'}[$vals[0]][$$outInfor{$outfile}{'dynamicMod'}{'found'}[$vals[0]]]=$j-$$outInfor{$outfile}{'dynamicMod'}{'found'}[$vals[0]]+1;
							$$outInfor{$outfile}{'dynamicMod'}{'mass'}[$vals[0]][$$outInfor{$outfile}{'dynamicMod'}{'found'}[$vals[0]]]=$$paraHash{'modification'}{'dynamic'}{'mass'}[$k];
						
							#print "\n$$outInfor{$outfile}{'dynamicMod'}{'found'}[$vals[0]]\:$i,$j,$k";
							last;
						}
					}
				}
				$j++;
			}
			
			$peptide =~ s/\W//g;
		}
		
		if (defined($$paraHash{'modification'}{'static'}{'count'}))		#static mod:
		{	#	print "\nstatic mod";	
			$$outInfor{$outfile}{'staticMod'}{'found'}[$vals[0]]=0;
			my @tmp=split(//,$peptide);
			for (my $j=0; $j<scalar(@tmp); $j++)
			{
				for (my $k=1; $k<=$$paraHash{'modification'}{'static'}{'count'}; $k++)
				{
					my $a=$$paraHash{'modification'}{'static'}{'AA'}[$k];
					if ($tmp[$j] eq $a)
					{
						$$outInfor{$outfile}{'staticMod'}{'found'}[$vals[0]]++;
						$$outInfor{$outfile}{'staticMod'}{'modIndex'}[$vals[0]][$$outInfor{$outfile}{'staticMod'}{'found'}[$vals[0]]]=$k;
						$$outInfor{$outfile}{'staticMod'}{'pos'}[$vals[0]][$$outInfor{$outfile}{'staticMod'}{'found'}[$vals[0]]]=$j+1;
						$$outInfor{$outfile}{'staticMod'}{'mass'}[$vals[0]][$$outInfor{$outfile}{'staticMod'}{'found'}[$vals[0]]]=$$paraHash{'modification'}{'static'}{'mass'}[$k];
						
						#print "\n$$outInfor{$outfile}{'staticMod'}{'found'}[$vals[0]]\:$i,$j,$k";
						last;
					}
				}
			}
			
		}
		
		$$outInfor{$outfile}{'peptide'}[$vals[0]]=$peptide;
		
		my @AA_KR=$peptide =~ /([KR])/g;
		my @AA_KRP=$peptide =~ /([KR]P)/g;
		$$outInfor{$outfile}{'num_tol_term'}[$vals[0]]=scalar(@AA_KR)-scalar(@AA_KRP);
			
		}
		else		#alternative proteins
		{
			$count++;
			$$outInfor{$outfile}{'altPro'}[$i][$count]=$vals[0];
			$$outInfor{$outfile}{'red'}[$i]=$count;
		}
	}
	
	if ($$outInfor{$outfile}{'hitShown'}==1) {$$outInfor{$outfile}{'deltacn'}[1]=1;}
	
	
	close(IN);
}

sub printPepXML
{
	my ($self,$seqPara, $outInfor, $run, $topHit, $outXML)=@_;
	my ($now_date,$now_time) = $self->gettime();
	my $Hydrogen_mass=1.007825032;	
	open(XML,">$outXML");

	print XML '<?xml version="1.0" encoding="Accoding to Out2XML(TPP v2.9 GALE rev.3, Build 200611091255)"?>',"\n";
	print XML '<msms_pipeline_analysis date="',"$now_date $now_time",'" ','summary_xml="',$run,$outXML,'">',"\n";
	print XML '<msms_run_summary base_name="',$run,'" raw_data_type="raw" raw_data="">',"\n";
	
	print XML '<sample_enzyme name="',$$seqPara{'enzyme'},'">',"\n";
	print XML '<specificity cut="',$self->enzymeStat($$seqPara{'enzyme'},'cut');
	unless ($self->enzymeStat($$seqPara{'enzyme'},'nocuts') eq '')
	{
		print XML '" no_cut="',$self->enzymeStat($$seqPara{'enzyme'},'nocuts');
	}
	print XML '" sense="',$self->enzymeStat($$seqPara{'enzyme'},'sense'),'"/>',"\n";
	print XML '</sample_enzyme>',"\n";
	
	print XML '<search_summary base_name="',$run,'" search_engine="SEQUEST" precursor_mass_type="monoisotopic" fragment_mass_type="monoisotopic" out_data_type="out" out_data=".out" search_id="">',"\n";
	print XML '<search_database local_path="',$$seqPara{'first_database_name'},'" type="AA"/>',"\n";
	if (defined($$seqPara{'modification'}{'dynamic'}{'count'}))
	{
		for (my $i=1; $i<=$$seqPara{'modification'}{'dynamic'}{'count'}; $i++)
		{
			print XML '<aminoacid_modification aminoacid="',$$seqPara{'modification'}{'dynamic'}{'AA'}[$i],'" massdiff="',$$seqPara{'modification'}{'dynamic'}{'massdiff'}[$i],'" mass="',$$seqPara{'modification'}{'dynamic'}{'mass'}[$i],'" variable="Y" symbol="',$$seqPara{'modification'}{'dynamic'}{'symble'}[$i],'"/>',"\n";
		}
	}
	if (defined($$seqPara{'modification'}{'static'}{'count'}))
	{
		for (my $i=1; $i<=$$seqPara{'modification'}{'static'}{'count'}; $i++)
		{
			print XML '<aminoacid_modification aminoacid="',$$seqPara{'modification'}{'static'}{'AA'}[$i],'" massdiff="',$$seqPara{'modification'}{'static'}{'massdiff'}[$i],'" mass="',$$seqPara{'modification'}{'static'}{'mass'}[$i],'" variable="N"/>',"\n";
		}
	}
	if (defined($$seqPara{'modification'}{'Nterm'}{'mass'}))
	{
		print XML '<terminal_modification terminus="n" massdiff="',$$seqPara{'modification'}{'Nterm'}{'massdiff'},'" mass="',$$seqPara{'modification'}{'Nterm'}{'mass'},'" variable="N" protein_terminus="N"/>',"\n";
	}
	if (defined($$seqPara{'modification'}{'Cterm'}{'mass'}))
	{
		print XML '<terminal_modification terminus="c" massdiff="',$$seqPara{'modification'}{'Cterm'}{'massdiff'},'" mass="',$$seqPara{'modification'}{'Cterm'}{'mass'},'" variable="N" protein_terminus="C"/>',"\n";
	}
	
	#OPTIONAL PARAMETERS
	foreach my $tmpTitle (@{$seqPara->{'optionalParams'}})
	{#print "$tmpTitle,$$seqPara{$tmpTitle}\n";
		if (defined($$seqPara{$tmpTitle}))
		{
			print XML "\<parameter name=\"$tmpTitle\" value=\"$$seqPara{$tmpTitle}\"\/\>\n";
		}
	}

	print XML '</search_summary>',"\n";
		
	opendir(DIR,$run); my @out = grep {/\.spout\Z/} readdir(DIR);
	closedir(DIR);
	my $count=0;
	for my $outfile (@out)
	{
		$count++;
		$outfile =~ s/(.*?)(\.spout)/$1$2/;  
		print XML "\<spectrum_query spectrum=\"$1\" start_scan=\"$$outInfor{$outfile}{'scan'}\" end_scan=\"$$outInfor{$outfile}{'scan'}\" precursor_neutral_mass=\"",$$outInfor{$outfile}{'MHmass'}-$Hydrogen_mass,
		"\" assumed_charge=\"$$outInfor{$outfile}{'charge'}\" index=\"$count\"\>\n";
		print XML '<search_result>',"\n";
		
		#unless (defined($$outInfor{$outfile}{'hitShown'})) {die "\n$outfile,$$outInfor{$outfile}{'hitShown'}\n";}
		for (my $i=1; $i<=$$outInfor{$outfile}{'hitShown'} and $i<=$topHit; $i++ )
		{
			#unless (defined($$outInfor{$outfile}{'peptide'}[$i]) and defined($$outInfor{$outfile}{'prAA'}[$i]) and defined($$outInfor{$outfile}{'nxAA'}[$i])) {die "\n$outfile\n";}
			print XML "\<search_hit hit_rank=\"",$i, '" peptide="',$$outInfor{$outfile}{'peptide'}[$i],'" peptide_prev_aa="',$$outInfor{$outfile}{'prAA'}[$i],'" peptide_next_aa="',$$outInfor{$outfile}{'nxAA'}[$i],
			'" protein="',$$outInfor{$outfile}{'protein'}[$i],'" num_tot_proteins="',$$outInfor{$outfile}{'red'}[$i]+1,'" num_matched_ions="',$$outInfor{$outfile}{'matchIons'}[$i],'" tot_num_ions="',
			$$outInfor{$outfile}{'totalIons'}[$i],'" calc_neutral_pep_mass="',$$outInfor{$outfile}{'expmass'}[$i] - $Hydrogen_mass,'" massdiff="',$$outInfor{$outfile}{'MHmass'}-$$outInfor{$outfile}{'expmass'}[$i],
			'" num_tol_term="',$$outInfor{$outfile}{'num_tol_term'}[$i],'" num_missed_cleavages="',$$seqPara{'misClv'},'" is_rejected="',0,"\"\>\n";
			
			#alternative_protein
			for (my $j=1; $j<=$$outInfor{$outfile}{'red'}[$i]; $j++)
			{
				unless (defined($$outInfor{$outfile}{'altPro'}[$i][$j])) {print "$outfile,red=$$outInfor{$outfile}{'red'}[$i],i=$i,j=$j\n";next;}
				print XML '<alternative_protein protein="',$$outInfor{$outfile}{'altPro'}[$i][$j],'"/>',"\n";
			}
			
			#modification
			if (defined($$seqPara{'modification'}{'Nterm'}{'mass'}) or defined($$seqPara{'modification'}{'Cterm'}{'mass'}) or 
			(defined($$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]) and $$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]) or 
			(defined($$outInfor{$outfile}{'staticMod'}{'found'}[$i]) and $$outInfor{$outfile}{'staticMod'}{'found'}[$i]) )
			{
				print XML '<modification_info';
				if (defined($$seqPara{'modification'}{'Nterm'}{'mass'}))
				{
					print XML  ' mod_nterm_mass="',$$seqPara{'modification'}{'Nterm'}{'mass'},"\"";
				}
				elsif (defined($$seqPara{'modification'}{'Cterm'}{'mass'}))
				{
					print XML  ' mod_cterm_mass="',$$seqPara{'modification'}{'Cterm'}{'mass'},"\"";
				}
				print XML "\>\n";
				
				if (defined($$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]))
				{#print "\nmark,$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]";
					for (my $j=1; $j<=$$outInfor{$outfile}{'dynamicMod'}{'found'}[$i]; $j++)
					{
						print XML '<mod_aminoacid_mass position="',$$outInfor{$outfile}{'dynamicMod'}{'pos'}[$i][$j],'" mass="',$$outInfor{$outfile}{'dynamicMod'}{'mass'}[$i][$j],'"/>',"\n";
					}
				}
				
				if (defined($$outInfor{$outfile}{'staticMod'}{'found'}[$i]))
				{#print "\nmark,$$outInfor{$outfile}{'staticMod'}{'found'}[$i]";
					for (my $j=1; $j<=$$outInfor{$outfile}{'staticMod'}{'found'}[$i]; $j++)
					{
						print XML '<mod_aminoacid_mass position="',$$outInfor{$outfile}{'staticMod'}{'pos'}[$i][$j],'" mass="',$$outInfor{$outfile}{'staticMod'}{'mass'}[$i][$j],'"/>',"\n";
					}
				}
				
				
				
				print XML '</modification_info>',"\n";
			}
			
			foreach my $tmpTitle (@{$seqPara->{'hitOptionalParams'}})
			{
				if (defined($$outInfor{$outfile}{$tmpTitle}[$i]))
				{
					print XML "\<search_score name=\"$tmpTitle\" value=\"$$outInfor{$outfile}{$tmpTitle}[$i]\"\/\>\n";
				}
			}
				
			print XML '</search_hit>',"\n";
		}
		
		print XML '</search_result>',"\n";
		print XML '</spectrum_query>',"\n";
	}
	
	print XML '</msms_run_summary>',"\n";
	print XML '</msms_pipeline_analysis>',"\n";
	
	close(XML);

}

1;

