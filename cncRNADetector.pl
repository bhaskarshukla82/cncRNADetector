use Getopt::Long;
#--------------------------configuration of external programs-------------------------------------------------
$cmd="";#could be modified according to users' computational environment
#----------------------------main interface that the program starts from--------------------------------------
&main(); 

#-----------------------------------subroutine begins----------------------------------------------------------
sub main(){
print"Welcome to use cncRNADetector 1.0: A Bioinformatics Pipeline for Coding and Non-Coding RNA Identification!\n";
name();
	my $USAGE = qq(
USAGE:
   lncRNADetector.pl -i <transcript.fasta> -p <protein.fasta> -k <housekeeping.fasta> -s <smallRNA.fasta> -o 
<output prefix> [-t <# of thread>] [-r <minimum lncRNA length>] [-f <maximum ORF length>] [-e <E-value of alignment>]
Options:
   -m <full/assembly/lncRNA/miRNA> select mode of pipeline 
   -p <pair/singular>  if use assembly workflow || default=pair
   -i   <transcript.fasta>
   -o <output prefix>
   -h help
   -t <int> number of thread for the computation || default=4
   -r <int> minimum lncRNA length || default=200
   -f <int> maximum potential ORF length of lncRNAs || default=120
   -m <int> number of mismatch in the alignment with smallRNA || default=0
   -e E-value of the alignment against protein database || default=1.0e-5
   
Requirement:
   Need to install ncbi_blast standalone package blastx, 2.2.28+,  Active Perl 5.28 and Cygwin for Windows platform in your local environment correctly and set the running path of external programs at the begining of the pipeline
   
Citation:
   Li L, Eichten SR, Shimizu R, Petsch K, Yeh CT, Wu W, Scanlon MJ, Yu JM, Schnable PS, Timmermans MCP, Springer NM, and Muehlbauer GJ: Genome-wide discovery and characterization of maize long non-coding RNAs (lncRNAs). Genome Biology, 2013, revised.
	);
	my ($inputfile,$hkfile,$srnafile,$profile,$outfile); #required parameters
	my ($mintslength,$maxorflength,$evalue,$thread);#optional parameters
	#
	my $parastr = &GetOptions(
							  "m=s{1}"=>\$mode,
							  "p=s{1}"=>\$dataset,
							  "i=s{1}"=>\$inputfile,
							  "k=s{1}"=>\$hkfile,
							  "s=s{1}"=>\$srnafile,
							  "q=s{1}"=>\$profile,
							  "o=s{1}"=>\$outfile,
							  "t=i{0,1}"=>\$thread,
							  "r=i{0,1}"=>\$mintslength,
							  "f=i{0,1}"=>\$maxorflength,
							  "e=s{0,1}"=>\$evalue
	);
	GetOptions ("help|?");

	unless ($parastr && defined($mode)) {
	   print $USAGE."\n";
	   exit 1;
	}else{
		print "----------------------------------------------------------\n\n";
		unless(-e $inputfile){
			print "input transcript sequence file-$inputfile does not exist\n";
			die("Pipeline halted\n");
		}
		
		#check the optional parameters
		unless(defined($thread)){
			$thread=4;
		}
		unless(defined($mintslength)){
			$mintslength=200;
		}
		unless(defined($maxorflength)){
			$maxorflength=120;
		}
		unless(defined($evalue)){
			$evalue=0.00001;
		}
		unless(defined($dataset)){
			$dataset="singular";
		}
		if ($mode eq "lncRNA"){

		print "----------------------------------------------------------\n\n";
		print "cncRNADetector 1.0 starts to work with lncRNA mode with the following parameters:\n";
		print "input transcript sequence file: $inputfile\n";
		print "number of thread for computation:\t$thread\n";
		print "minimum lncRNA length:\t$mintslength\n";
		print "maximum potential ORF length of lncRNAs:\t$maxorflength\n";
		print "E-value of the alignment against protein database:\t$evalue\n";
		print "----------------------------------------------------------\n\n";
	
#my $inputfile;
#my $parastr = &GetOptions("i=s{1}"=>\$inputfile);
my @folder=split(".fa",$inputfile);
#print "@folder[0]";
my $folder=@folder[0];
$file=$folder;
#@file=split("\/",$file);
#$file=@file[1];
swissblast($inputfile,$evalue,$folder);
pataa($inputfile,$evalue,$folder);
codingpotiential($folder);
orf_length($mintslength,$maxorflength,$folder);

		}
		elsif ($mode eq "assembly"){
			print "----------------------------------------------------------\n\n";
		print "cncRNADetector 1.0 starts to work with Assembly mode with the following parameters:\n";
		print "input transcript sequence file: $inputfile\n";
			my @file=split(".SRA",$inputfile);
#print "@folder[0]";
my $file=@file[0];
@filename=split /\\/,$file;
$filenames=pop(@filename);
$cmd="fastq-dump --split-files  ".$filenames;
system($cmd);
if($dataset eq "singular"){
$cmd ="perl trim_galore -q 30 --gzip --length 35 --paired ".$filenames."_1.fastq  ".$filenames."_2.fastq";
system($cmd);
$cmd="/home/shivam/Downloads/jakeesir/trinityrnaseq-v2.14.0.FULL/trinityrnaseq-v2.14.0/Trinity --seqType fq  --single /home/shivam/TrimGalore-master//".$filenames."_1_trimmed.fq.gz --max_memory 7G  --CPU 6  --full_cleanup --output /home/shivam/Desktop/cncDector/trinity";
system($cmd);
}


}
else {
	print "select correct option in mode parameter \n";
}
	}
#python tool/CPC2/bin/CPC2.py  -i result/Datura_innoxia_Output_lenfilterleft.fasta -o result/Datura_innoxia_CPC2.txt

#python tool/PLEK/PLEK.py  -fasta result/Datura_innoxia_Output_lenfilterleft.fasta -out result/Datura_innoxia_PLEK.txt

#python2 tool/lgc/LGC-1.0.py  result/Datura_innoxia_Output_lenfilterleft.fasta  result/Datura_innoxia_LGC.txt

#python3 ./src/mipepid.py ./demo_files/sample_seqs.fasta ./demo_files/MiPepid_results_on_sample_seqs.csv
}
sub swissblast(){
	my($inputfile,$evalue,$folder)=@_;
print "\n\n--------------------------------swiss blastx run---------------------------------\n\n";
my $cmd="blastx -query $inputfile  -db db/swiss/swissport  -evalue $evalue -out ".$folder."_blast_swissoutput.txt -qcov_hsp_perc 90 -max_hsps 1";
system($cmd);
open my $keywords,    '<','search.txt' or die "Can't open keywords: $!";
open my $search_file, "<  ".$folder."_blast_swissoutput.txt"   or die "Can't open search file: $!";
my $keyword_or = join '|', map {chomp;qr/\Q$_\E/} <$keywords>;
my $regex = qr|\b($keyword_or)\b|;
my @NO_HIT;
while (<$search_file>)
{
    while (/$regex/g)
    {
       		push @NO_HIT, $.;
			#print "$.: $1\n";
			
			
    }
}
#print "@NO_HIT\n";
my $t= scalar(@NO_HIT);
#print "$t\n"  total not match;
open(hh,">",$folder."total_no_seq_not_match__swiss.txt");
print hh " total  no hit sequences in swiss : $t \n";

close(hh);
open my $keywords2,    '<', 'search2.txt' or die "Can't open keywords: $!";
open my $search_file2, '<', $folder.'_blast_swissoutput.txt'   or die "Can't open search file: $!";
my $keyword_or2 = join '|', map {chomp;qr/\Q$_\E/} <$keywords2>;
my $regex2 = qr|\b($keyword_or2)\b|;
my @QUERY;
while (<$search_file2>)
{
    while (/$regex2/g)
    {
       		push @QUERY, $.;
			#print "$.: $1\n";
			
    }
}
#print "@QUERY\n";
my @psn;
for my $i (0..$t-1){
my @suny2=grep{$_ <$NO_HIT[$i]}@QUERY;
my $sun3=pop(@suny2);
push @psn, $sun3; 
}
#print "@psn\n";

#Finding out the leftsequence after BLASTX where No hits found
my $size=scalar(@psn)-1;
my @leftsequence;
open (FILE, "<", $folder.'_blast_swissoutput.txt') or die "could not open filename";
my @lines = <FILE>;
close FILE;
while($size>=0)
{
#print $lines[@psn[$size]-1];
#push @leftsequence,$lines[@psn[$size]-1];
push @leftsequence,substr $lines[$psn[$size]-1],7;
$size--;
}
#print "@leftsequence\n";
foreach $seq(@leftsequence){
@leftsequence_acc=split(" ",$seq);
push(@leftsequence_acc_array,">@leftsequence_acc[0]");
}
#print @leftsequence_acc_array;
$inputfile=$inputfile;
open(r,"<", $inputfile);
@a1=<r>;
$fasta=join("",@a1);
@fasta=split(">",$fasta);
foreach $file(@fasta)
{
@fasta_lines=split("\n",$file);
#print " @fasta_lines[0] \n";
$fasta_first_line=join("",@fasta_lines[0]);
@fasta_first_line=split(" ",$fasta_first_line);
shift(@fasta_lines);
$seq=join("",@fasta_lines);
#print " $seq  \n";
push(@seq, $seq); 
#print "@fasta_first_line[0] \n";
push(@acc_id,">@fasta_first_line[0]");
}
shift @acc_id;
shift @seq;
# make diect
$j=0;
for($i=0;$i<=scalar @acc_id ;$i++)
{
push(@dic,@acc_id[$i]) ;
push(@dic,@seq[$i])
}
%hash1=@dic;


foreach $key(keys %hash1)
{
foreach $index(@leftsequence_acc_array){
if($key eq $index)
{
push(@not_match,$index);
}
}}
#write sequence not found in pattaa
open(w,"> ".$folder."_not_in_swiss.txt");
foreach $a(@not_match)
{
$value=$hash1{$a};
print w "$a \n";
@$value=split("",$value);

$value=~s/^\s+|\s+$//g;
print w"$value \n";

}
close(w);
print "\n\n--------------------------------swiss blastx  complete---------------------------------\n\n";
}
sub pataa(){
	my($inputfile,$evalue,$folder)=@_;
print "\n\n--------------------------------pataa blastx  start---------------------------------\n\n";
$inputfile=$folder."_not_in_swiss.txt";
my $cmd="blastx -query $inputfile  -db db/pataa/pataa  -evalue $evalue -out ".$folder."_blast_pataaoutput.txt -qcov_hsp_perc 90 -max_hsps 1";
system($cmd);
open my $keywords2,    '<', 'search.txt' or die "Can't open keywords: $!";
open my $search_file2, '<', $folder.'_blast_pataaoutput.txt'   or die "Can't open search file: $!";
my $keyword_or2 = join '|', map {chomp;qr/\Q$_\E/} <$keywords2>;
my $regex2 = qr|\b($keyword_or2)\b|;
my @NO_HIT2;
while (<$search_file2>)
{
    while (/$regex2/g)
    {
       		push @NO_HIT2, $.;
			#print "$.: $1\n";
			
			
    }
}
#print "@NO_HIT\n";
my $t2= scalar(@NO_HIT2);
#print "$t\n"  total not match;
open(hh,">",$folder."total_no_seq_not_match__pataa.txt");
print hh " total  no hit sequences in pataa : $t2 \n";

close(hh);
open my $keywords2,    '<', 'search2.txt' or die "Can't open keywords: $!";
open my $search_file2, '<', $folder.'_blast_pataaoutput.txt'   or die "Can't open search file: $!";
my $keyword_or2 = join '|', map {chomp;qr/\Q$_\E/} <$keywords2>;
my $regex2 = qr|\b($keyword_or2)\b|;
my @QUERY2;
while (<$search_file2>)
{
    while (/$regex2/g)
    {
       		push @QUERY2, $.;
			#print "$.: $1\n";
			
    }
}
#print "@QUERY\n";
my @psn2;
for my $i (0..$t2-1){
my @suny2=grep{$_ <$NO_HIT2[$i]}@QUERY2;
my $sun3=pop(@suny2);
push @psn2, $sun3; 
}
#print "@psn\n";

#Finding out the leftsequence after BLASTX where No hits found
my $size=scalar(@psn2)-1;
my @leftsequence2;

open (FILE, "<", $folder.'_blast_pataaoutput.txt') or die "could not open filename";
my @lines = <FILE>;
close FILE;
while($size>=0)
{
#print $lines[@psn[$size]-1];
#push @leftsequence,$lines[@psn[$size]-1];
push @leftsequence2,substr $lines[$psn2[$size]-1],7;
$size--;
}
#print "@leftsequence\n";
foreach $seq(@leftsequence2){
@leftsequence_acc=split(" ",$seq);
push(@leftsequence_acc_array2,">@leftsequence_acc[0]");
}
#print @leftsequence_acc_array;
$inputfile2=$folder.'_not_in_swiss.txt';
open(r,"<", $inputfile2);
@a1=<r>;
$fasta=join("",@a1);
@fasta=split(">",$fasta);
foreach $file(@fasta)
{
@fasta_lines=split("\n",$file);
#print " @fasta_lines[0] \n";
$fasta_first_line=join("",@fasta_lines[0]);
@fasta_first_line=split(" ",$fasta_first_line);
shift(@fasta_lines);
$seq=join("",@fasta_lines);
#print " $seq  \n";
push(@seq2, $seq); 
#print "@fasta_first_line[0] \n";

push(@acc_id2,">".@fasta_first_line[0]);
}
shift @acc_id2;
shift @seq2;
# make diect
$j=0;
for($i=0;$i<=scalar @acc_id2 ;$i++)
{
push(@dic2,@acc_id2[$i]) ;
push(@dic2,@seq2[$i])
}
%hash2=@dic2;
foreach $key(keys %hash2)
{
foreach $index(@leftsequence_acc_array2){
if($key eq $index)
{
push(@not_match2,$index);
}
}}
#write sequence not found in pattaa
open(w,"> ".$folder."_not_in_pataa.fasta");
foreach $a(@not_match2)
{
$value=$hash2{$a};
$a=~ s/^\s+|\s+$//g;
chomp($a);
print w "$a\n";
@value=split("\n",$value);
$value=join(" ",@value);
print w "$value\n";
}

close(w);
#write coding-sequence  found in pattaa
@id=keys %hash1;
foreach $id(@id)
{
push(@patta_match2, $id) unless grep { $id eq $_ } @not_match2;
}	
open(ww,"> ".$folder."coding.txt");
foreach $a(@patta_match2)
{
$value=$hash1{$a};
print ww "$a \n";
print ww"$value \n";
}
close(ww);
#print "\n Total no of  seq found in pataa db :  ";
##print "\n\n";
#print "\n Total no of  seq found in pataa db :  ";
#print scalar @patta_match;
#print "\n\n";
print "\n\n--------------------------------pataa blast complete---------------------------------\n\n";
}
sub orf_length(){
my ($mintslength,$maxorflength,$folder)=@_;
$files=$folder."_noncoding.fasta";
print $files;
my  %oriseqhash=getOriSeqinfo($files);
	print "processing transcript length/ORF filter...\n";
	print "----------------------------------------------------------\n\n";
	#Step 1 & 2 Code Starts: (Transcript Size Selection)/ORF filter
	my $tother=$folder."_other_ncRNA.txt";
	my $tslorf=$folder."_RNAlen_ORFlen_info.txt";
	my $lenfilterleftseqfile=$folder."_lncRNA.fasta";
        open otherout,">$tother"|| die "Cannot create result file:$!";
	open LENOUT,">$tslorf" || die "Cannot create result file:$!";
	open LEFTOUT,">$lenfilterleftseqfile" || die "Cannot create result file:$!";
	print LENOUT "seqname\ttranscriptlen\tORF1\tORF2\tORF3\trevORF1\trevORF2\trevORF3\n";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		my $toutstr=$tkey;
		$toutstr.="\t".length($tvalue);
		my $checkstr=checkprotein(translate_frame($tvalue, 1));
		my $checkstr2=checkprotein(translate_frame($tvalue, 2));
		my $checkstr3=checkprotein(translate_frame($tvalue, 3));
		$toutstr.="\t".length($checkstr)."\t".length($checkstr2)."\t".length($checkstr3);
		my $revcom = revcom($tvalue);
		my $checkstr4=checkprotein(translate_frame($revcom, 1));
		my $checkstr5=checkprotein(translate_frame($revcom, 2));
		my $checkstr6=checkprotein(translate_frame($revcom, 3));
		$toutstr.="\t".length($checkstr4)."\t".length($checkstr5)."\t".length($checkstr6);
		print LENOUT $toutstr."\n";
		if(length($checkstr)>$maxorflength || length($checkstr2)>$maxorflength || length($checkstr3)>$maxorflength || length($checkstr4)>$maxorflength || length($checkstr5)>$maxorflength || length($checkstr6)>$maxorflength || length($tvalue)<$mintslength){
#print LENOUT ">".$tkey."\n".$tvalue."\n"; write non coding other
         print otherout ">".$tkey."\n".$tvalue."\n";
			#$oriseqhash{$tkey}="";
		}else{
#print $tkey."\n".$tvalue."\n";
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	close(LENOUT);
#close(otherout);
#Step 1 & 2 Closed : (Transcript Size Selection)/ORF filter	

	print "transcript length/ORF filter completed\n";
	print "----------------------------------------------------------\n\n";
}
sub codingpotiential(){
	my($folder)=@_;
$cmd="python tool/CPC2/bin/CPC2.py  -i ".$folder."_not_in_pataa.fasta -o ".$folder."_CPC2";
system($cmd);
$cmd="python tool/PLEK/PLEK.py  -fasta ".$folder."_not_in_pataa.fasta -out ".$folder."_innoxia_PLEK.txt";
system($cmd);
$cmd="python2 tool/lgc/LGC-1.0.py  ".$folder."_not_in_pataa.fasta  ".$folder."_LGC.txt";
system($cmd);
open(hhh,"<".$folder."__not_in_pataa.fasta");
@filedata1=<hhh>;
%hash11=@filedata1;
#print"enter file name  pleck:";
$plek=$folder."_innoxia_PLEK.txt";
chomp($plek);
open(FH, "< $plek") or die"Can't open file: $!";
while(<FH>)
{
@lines=split("\t",$_);
$coding=@lines[0]."\t".@lines[1];
$id=@lines[2];
$id=~ s/^\s+|\s+$//g;
push(@array,$id);
push(@array,$coding);
}
%plek=@array;
#print"enter file name  cpc :";
$cpc=$folder."_CPC2.txt";
chomp($cpc);
open(FH, "< $cpc") or die"Can't open file: $!";
@file=<FH>;
shift(@file);
foreach $b(@file)
{
@lines=split("\t",$b);
$id=@lines[0];
$id=~ s/^\s+|\s+$//g;
$coding=@lines[6]."\t".@lines[7];
push(@cpcarray,$id);
push(@cpcarray,$coding);
}
%cpc=@cpcarray;
#print"enter file name  lgc:";
$lgc=$folder."_LGC.txt";
chomp($lgc);
open(FH, "< $lgc") or die"Can't open file: $!";
@file=<FH>;
for($i=0;$i<11;$i++)
{
	shift(@file);
}
foreach $a(@file)
{
@lines=split("\t",$a);
$id=@lines[0];
$id=~ s/^\s+|\s+$//g;
$coding=@lines[3]."\t".@lines[4];
push(@lgcarray,$id);
push(@lgcarray,$coding);
}
%lgc=@lgcarray;
foreach $key( keys %cpc)
{
	$key2=">".$key;
	$plek_value=$plek{$key2};
	@plek_noncoding=split("\t",$plek_value);
	$plek_P=@plek_noncoding[1];
	push(@plek_P,$key2);
	push(@plek_P,$plek_P);
	$plek_noncoding=@plek_noncoding[0];
	$plek_noncoding=~ s/^\s+|\s+$//g;
	$cpc_value=$cpc{$key};
	@cpc_noncoding=split("\t",$cpc_value);
	$cpc2_P=@cpc_noncoding[0];
	push(@cpc2_P,$key2);
	push(@cpc2_P,$cpc2_P);
	$cpc_noncoding=@cpc_noncoding[1];
	$cpc_noncoding=~ s/^\s+|\s+$//g;
        $plek_noncoding=lc($plek_noncoding);
	@plek_noncoding=split("-",$plek_noncoding);
	$plek_noncoding=join("",@plek_noncoding);
	$lgc_value=$lgc{$key};
	@lgc_noncoding=split("\t",$lgc_value);
	$lgc_P=@lgc_noncoding[0];
	push(@lgc_P,$key2);
	push(@lgc_P,$lgc_P);
	$lgc_noncoding=@lgc_noncoding[1];
	$lgc_noncoding=~ s/^\s+|\s+$//g;
        $lgc_noncoding=lc($lgc_noncoding);
	@lgc_noncoding=split("-",$lgc_noncoding);
	$lgc_noncoding=join("",@lgc_noncoding);
	
	if((($plek_noncoding eq "noncoding") && ($cpc_noncoding eq "noncoding")) || (($cpc_noncoding eq "noncoding") && ($lgc_noncoding eq "noncoding")) || (($plek_noncoding eq "noncoding") && ($lgc_noncoding eq "noncoding")))
	{#print "$plek_noncoding  \t $cpc_noncoding  \t  $lgc_noncoding  \n";
		push(@common,$key)
	}
	else
	{
		push(@not_common,$key);
		#print  " not match   :$key \n";
	}
}
%plek_P=@plek_P;
%cpc2_P=@cpc2_P;
%lgc_P=@lgc_P;
	print "-------------------PLEK/CPC/LGC completed----------------\n";
	print "----------------------------------------------------------\n\n";

open(FHH, "> ".$folder."_noncoding.fasta") or die"Can't open file: $!";
foreach $lnc(@common){
$key=">".$lnc;
$value=$hash1{$key};
$v1=$plek_P{$key};
$v2=$cpc2_P{$key};
$v3=$lgc_P{$key};
print FHH ">$lnc \t CPC2:$v2 \t LGC:$v3 \t PLEK:$v1 \n";
print FHH "$value \n";
}
print FHH ">\n";
open(ww,">>".$folder."coding.txt");
print ww "coding potentail \n";
foreach $lnc(@not_common){
$key=">".$lnc;
$value=$hash1{$key};
print ww ">$lnc  \n";
print ww "$value \n";
}
}
sub getOriSeqinfo(){
#print " enter file name ";
    #my $namefiles=$file."_not_in_pataa.txt";
	my $orifile=shift;
	
	my %resarr;
	
	open ORI,$orifile or die "Cannot open $orifile:$!";
	my $tstr;
	while($tstr=<ORI>){
		if(trim($tstr) ne "" && $tstr=~/>/){
			my @tarr=split(/\s+/,trim($tstr));
			my $tseqstr="";
			my $nstr="";
			while ($nstr=<ORI>){
				if($nstr=~/>/){
					my $tname=substr(trim($tarr[0]),1);
					$resarr{$tname}=$tseqstr;
					seek(ORI,-length($nstr),1);
					last;
				}else{
					$tseqstr.=trim($nstr);
				}
			}
		}
	}
	close(ORI);

	return %resarr;
}
sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }#else{

#            print STDERR "Bad codon \"$codon\"!!\n";
#            exit;
#    }
}

sub dna2peptide {

    my($dna) = @_;

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
	
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}

sub translate_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    # To make the subroutine easier to use, you won't need to specify
    #  the end point-it will just go to the end of the sequence
    #  by default.
    unless($end) {
        $end = length($seq);
    }

    # Finally, calculate and return the translation
    return dna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}

# revcom 
#
# A subroutine to compute the reverse complement of DNA sequence

sub revcom {

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}

sub checkprotein(){
	my $tstr=shift;
	my $resstr="";
	my @tarr=split(/\*/,trim($tstr));
	#Code modified here 
	
	my $larr=rindex(trim($tstr),"*");
	push (@tarr,  substr trim($tstr), $larr + 1, length(trim($tstr))-$larr);
	# ................
	
	my $maxlen=0;
	my $maxstr="NA";
	my $tpro="";
	pop(@tarr);
	foreach $tpro(@tarr){
		if(trim($tpro) ne ""){
			if(index($tpro,"M")!=-1){
				if($maxlen<length(trim($tpro))-index(trim($tpro),"M")){
					$maxlen=length(trim($tpro))-index(trim($tpro),"M");
					$maxstr=substr(trim($tpro),index(trim($tpro),"M"));
				}
			}
		}
	}
	$resstr=$maxstr;
	return $resstr;
}

sub trim(){
	my $string=shift;
	$string=~s/^\s+//;
	$string=~s/\s+$//;
	return $string;}
sub name(){

print" ----------------------------------------------------------------------------------------------- \n";
print"|                                                                                               |\n";
print"|                    ####   #    #  ####   ####            #                 #                  |\n";
print"|                    #   #  ##   # #    #  #   #           #                 #                  |\n";
print"|  #### #   #  ##### #    # ##   # #    #  #    #  ####  ##### ####   #### #####  ###   #  ###  |\n";
print"| #     #   # #      #####  # #  # #    #  #     # #   #   #   #   # #       #   #   #  # #   # |\n";
print"| #     ##  # #      ##     #  # # ######  #    #  ####    #   ####  #       #   #   #  ##      |\n";
print"| #     #  ## #      # #    #   ## #    #  #   #   #       #   #     #       #   #   #  #       |\n";
print"|  #### #   #  ##### #  ### #    # #    #  ####    ####    ### ####   ####   ###  ###   #       |\n";
print"|                                                                                               |\n";
print" -----------------------------------------------------------------------------------------------\n";
}
