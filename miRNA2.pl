use Getopt::Long;
my $inputfile;
my $parastr = &GetOptions("i=s{1}"=>\$inputfile);
my @folder=split(".fa",$inputfile);
my @folder2=split /\\/,$inputfile;
$folder_name=@folder2[0];
#print "@folder[0]";
my $folder=@folder[0];
#print "$folder";
#shift(@folder,$folder);
$path="C:\\xampp\\htdocs\\lncRNA\\".$folder_name;
my $path2="C:\\xampp\\htdocs\\lncRNA\\".$folder;
print "\n\n--------------------------------tRNA run---------------------------------\n\n";
$db="C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\db\\tRNA\\tRNA";
my $cmd="C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\blastn -db ".$db." -query C:\\xampp\\htdocs\\lncRNA\\".$folder."_other_ncRNA.fasta -evalue 0.00001 -out C:\\xampp\\htdocs\\lncRNA\\".$folder_name."\\miRNA\\_blast_tRNAoutput.txt";
system($cmd);
#print "\n".$cmd."\n";
open my $keywords, '<','search.txt' or die "Can't open keywords: $!";
open my $search_file, "< ".$path."\\miRNA\\_blast_tRNAoutput.txt" or die "Can't open search file: $!";
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
close(hh);
open my $keywords2, '<', 'search2.txt' or die "Can't open keywords: $!";
open my $search_file2, "< ".$path."\\miRNA\\_blast_tRNAoutput.txt" or die "Can't open search file: $!";
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
open (FILE, "< ".$path."\\miRNA\\_blast_tRNAoutput.txt") or die "could not open filename";
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
$inputfile=$path2."_other_ncRNA.fasta";
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
#print " $seq  \n";
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
open(w,"> ".$path."\\miRNA\\_not_in_tRNA.txt");
foreach $a(@not_match)
{
$value=$hash1{$a};
print w "$a \n";
print w"$value \n";
}
close(w);

print "\n\n--------------------------------rRNA blastn start---------------------------------\n\n";
$inputfile=$path."\\miRNA\\_not_in_tRNA.txt";
my $cmd="C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\blastn -db C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\db\\tRNA\\rRNA -query $inputfile -evalue 0.00001 -out ".$path."\\miRNA\\_blast_rRNAoutput.txt";
system($cmd);
open my $keywords2, '<', 'search.txt' or die "Can't open keywords: $!";
open my $search_file2, '<', $path."\\miRNA\\_blast_rRNAoutput.txt" or die "Can't open search file: $!";
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
#print "$t\n"  total not match;
close(hh);
open my $keywords2, '<', 'search2.txt' or die "Can't open keywords: $!";
open my $search_file2, '<', $path."\\miRNA\\_blast_rRNAoutput.txt"or die "Can't open search file: $!";
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

open (FILE, "<", $path."\\miRNA\\_blast_rRNAoutput.txt") or die "could not open filename";
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
$inputfile2=$path."\\miRNA\\_not_in_tRNA.txt";
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
#print " $seq  \n";
push(@seq2, $seq);
#print "@fasta_first_line[0] \n";
push(@acc_id2,">@fasta_first_line[0]");
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
open(w,"> ".$path."\\miRNA\\_not_in_rRNA.txt");
foreach $a(@not_match2)
{
$value=$hash2{$a};
print w "$a \n";
print w"$value \n";
}

close(w);
#write coding-sequence  found in pattaa
@id2=keys %hash2;
foreach $id(@id2)
{
push(@patta_match2, $id) unless grep { $id eq $_ } @not_match2;
}
open(ww,"> ".$path."\\miRNA\\t-rRNA.txt");
foreach $a(@patta_match2)
{
$value=$hash2{$a};
print ww "$a \n";
print ww"$value \n";
}
close(ww);
print "\n\n--------------------------------miRNA  blastn  start---------------------------------\n\n";
$inputfile=$path."\\miRNA\\_not_in_rRNA.txt";
my $cmd="C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\blastn -db C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\db\\mature\\miRNA -query ".$inputfile." -out ".$path."\\miRNA\\_blast_miRNAoutput.txt -evalue 0.1 -word_size 4   -outfmt ".'"0 mismatch 4"';
system($cmd);
open my $keywords2, '<', 'search.txt' or die "Can't open keywords: $!";
open my $search_file2, '<', $path."\\miRNA\\_blast_miRNAoutput.txt" or die "Can't open search file: $!";
my $keyword_or2 = join '|', map {chomp;qr/\Q$_\E/} <$keywords2>;
my $regex2 = qr|\b($keyword_or2)\b|;
my @NO_HIT3;
while (<$search_file2>)
{
while (/$regex2/g)
{
push @NO_HIT3, $.;
#print "$.: $1\n";
}
}
#print "@NO_HIT\n";
my $t3= scalar(@NO_HIT3);
#print "$t\n"  total not match;


close(hh);
open my $keywords2, '<', 'search2.txt' or die "Can't open keywords: $!";
open my $search_file2, '<', $path."\\miRNA\\_blast_miRNAoutput.txt" or die "Can't open search file: $!";
my $keyword_or2 = join '|', map {chomp;qr/\Q$_\E/} <$keywords2>;
my $regex2 = qr|\b($keyword_or2)\b|;
my @QUERY3;
while (<$search_file2>)
{
while (/$regex2/g)
{
push @QUERY3, $.;
#print "$.: $1\n";
}
}
#print "@QUERY\n";
my @psn3;
for my $i (0..$t3-1){
my @suny3=grep{$_ <$NO_HIT3[$i]}@QUERY3;
my $sun3=pop(@suny3);
push @psn3, $sun3;
}
#print "@psn\n";

#Finding out the leftsequence after BLASTX where No hits found
my $size=scalar(@psn3)-1;
my @leftsequence3;

open (FILE, "<", $path."\\miRNA\\_blast_miRNAoutput.txt") or die "could not open filename";
my @lines = <FILE>;
close FILE;
while($size>=0)
{
#print $lines[@psn[$size]-1];
#push @leftsequence,$lines[@psn[$size]-1];
push @leftsequence3,substr $lines[$psn3[$size]-1],7;
$size--;
}
#print "@leftsequence\n";
foreach $seq(@leftsequence3){
@leftsequence_acc=split(" ",$seq);
push(@leftsequence_acc_array3,">@leftsequence_acc[0]");
}
#print @leftsequence_acc_array;
$inputfile2=$path.'\\miRNA\\_not_in_rRNA.txt';
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
#print " $seq  \n";
push(@seq3, $seq);
#print "@fasta_first_line[0] \n";
push(@acc_id3,">@fasta_first_line[0]");
}
shift @acc_id3;
shift @seq3;
# make diect
$j=0;
for($i=0;$i<=scalar @acc_id3 ;$i++)
{
push(@dic3,@acc_id3[$i]) ;
push(@dic3,@seq3[$i])
}
%hash3=@dic3;
foreach $key(keys %hash3)
{
foreach $index(@leftsequence_acc_array3){
if($key eq $index)
{
push(@not_match3,$index);
}
}}
#write sequence not found in pattaa
open(w,"> ".$path."\\miRNA\\other_RNA.txt");
foreach $a(@not_match3)
{
$value=$hash3{$a};
print w "$a \n";
print w"$value \n";
}
close(w);
#write coding-sequence  found in pattaa
@id3=keys %hash3;
foreach $id(@id3)
{
push(@patta_match3, $id) unless grep { $id eq $_ } @not_match3;
}
open(ww,"> ".$path."\\miRNA\\_pri-miRNA.fasta");
foreach $a(@patta_match3)
{
$value=$hash3{$a};
if($value ne ''){
print ww "$a \n";
print ww"$value \n";
}
}
close(ww); 
if (0) {
    #print "-----This is an empty file----";-z $path."\\miRNA\\_pri-miRNA.fasta"
	
	open(hh,"> ".$path."\\miRNA\\miRNA.txt");
}
else{
print"---------RNAfold start--------\n";
$path2="C:\\xampp\\htdocs\\lncRNA\\".$folder;

$cmd="C:\\ViennaRNA\\RNAfold  < ".$path."\\miRNA\\_pri-miRNA.fasta > ".$path."\\miRNA\\_fold.out   --noPS ";

system($cmd);
print"\n---------RNAfold complete--------";
open(h,"< ".$path."\\miRNA\\_fold.out");
@fold=<h>;
$size=scalar @fold;
for ($i=0;$i<=$size;$i++)
{$fold_key=@fold[$i];
$fols_seq=@fold[$i+1].@fold[$i+2];
push(@fold_dataa,$fold_key);
#print $fold_key."\n";
push(@fold_dataa,$fols_seq);
$i+=2;
}
%fold_data=@fold_dataa;
foreach $fold_key(keys %fold_data)
{
$value=$fold_data{$fold_key};
@line=split("\n",$value);
$seq=@line[0];
$a= @line[1];
@word=split("",$a);
@pos=();
while($a=~m/\)\.*\(/g)
{
$pos=pos($a);
push(@pos,$pos);
#print $pos."\n";
}
$size=scalar(@word);
$j=0;
$k=0;
$last=@pos[0];
$id=0;
for($i=0;$i<=$size;$i++)
{
if($i == @pos[$j])
{
$temp=substr($seq,$k,$last);
$temp2=substr($a,$k,$last);
chomp($fold_key);
push(@str,$fold_key."-".$id);
push(@str,$temp);
$id++;
$k=@pos[$j]-2;
$j++;
$f=@pos[$j];
$last=$f-$k;
#print $id."\t\n";
#print $temp2;
}
#print @word[$i];
}
}
%hash=@str;
open(miRNA,"> ".$path."\\miRNA\\pre-miRNA.fasta");
foreach $key(keys %hash)
{
$value=$hash{$key};
if($value ne '')
{
print miRNA $key."\n";
print miRNA "$value \n";
}
else
{
@key=split("-",$key);
$key=@key[0];
$key=$key."\n";
$value=$fold_data{$key};
#print "\n".$key;
@line=split("\n",$value);
$seq=@line[0];
chomp($key);
if($seq ne ''){
print miRNA $key."-0"."\n";
print miRNA "$seq \n";
}}
}

print"-------validation stem-loop----------------";
print"---------RNAfold start--------\n";
#$cmd="C:\\ViennaRNA_Package\\RNAfold  --noPS < ".$path."\\miRNA\\_pri-miRNA.fasta > ".$path."\\miRNA\\_fold.out";
$cmd="C:\\ViennaRNA\\RNAfold  --noPS < ".$path."\\miRNA\\pre-miRNA.fasta > ".$path."\\miRNA\\_pre-miRNA_fold.out";
system($cmd);
print"\n---------RNAfold complete--------";
open(ff,"< $inputfile");
@fold=<ff>;
$size=scalar @fold;
for($i=0;$i<=$size;$i++)
{
$keyss=@fold[$i];
chomp($keyss);
$keyss=~s/^\s+|\s+$//g;
push (@foldd,$keyss);
push (@foldd,@fold[$i+1]);
$i=$i+1;
}
%fildataa=@foldd;
$filename=$path."\\miRNA\\_pre-miRNA_fold.out";

open(h,"< $filename");
@fold=<h>;
$size=scalar @fold;
for ($i=0;$i<=$size;$i++)
{
$key=@fold[$i+2];
@mef=split(" ",$key);
$mef=pop @mef;
@mef=split(/\(/,$mef);
$mef=pop @mef;
@mef=split(/\)/,$mef);
if (@mef[0]< -9){

push(@keyssmi,@fold[$i]);
push(@keyssmi,@fold[$i+1]);
push(@keyss,@fold[$i]);
push(@keyss,@fold[$i+1].@fold[$i+2]);
push(@keyssmef,@fold[$i]);
push(@keyssmef,@mef[0]);
@b=split("",@fold[$i+1]);
$total= scalar @b;
$count=0;
foreach $b(@b)
{ 
  if(($b eq 'G')||($b eq 'g')||($b eq 'c')||($b eq 'C'))
{
   $count++; 
}}
$gc= ($count/$total)*100;
push(@keyssgc,@fold[$i]);
push(@keyssgc,$gc);
}
$i=$i+2;
}
%precusor=@keyss;
%mefs=@keyssmef;
%gc=@keyssgc;
%pre_miRNa=@keyssmi;
$filepre=$path."\\miRNA\\_precusor.fasta";
open(h,"> $filepre");
foreach $keys(keys %precusor)
{
$value=$precusor{$keys};
$gc=$gc{$keys};
$mef=$mefs{$keys};
$value_seq=$pre_miRNa{$keys};
chomp($keys);
@idd=split("-",$keys);
$idd=@idd[0];
$idd=~s/^\s+|\s+$//g;
#$idd=$idd."\n";
$values=$fildataa{$idd};
$value_seq=~s/U/T/g;
chomp($value_seq);
#$index = index ($string, 'the');
$index = index($values, $value_seq);
$size=length($value);
$end=$index+$size;
$position=$index." ".$end;
#print $position;
print h "$keys\tgc:$gc\tmef:$mef\t$position\n";
print h "$value\n";
}
open(hh,"> miRNA/".$file."pre_miRNA.fasta");
foreach $keys(keys %pre_miRNa)
{
$value=$pre_miRNa{$keys};
print hh "$keys";
print hh "$value\n";
}
$inputfile=$path."\\miRNA\\pre-miRNA.fasta";
#my $cmd="C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\blastn -db C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\db\\mature\\miRNA -query ".$inputfile." -out ".$path."\\miRNA\\_blast_miRNAoutput.txt -evalue 0.1 -word_size 7   -outfmt ".'"0 mismatch 3"';
my $cmd="C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\blastn -db C:\\ncbi\\ncbi-blast-2.12.0+-x64-win64\\ncbi-blast-2.12.0+\\bin\\db\\mature\\miRNA -query ".$inputfile." -evalue 0.1 -word_size 7 -out $path\\miRNA\\_blast_Pre-miRNAoutput.txt -outfmt ".'"  0 mismatch 3"';
system($cmd);
open my $search_file2, '<', $path.'\\miRNA\\_blast_Pre-miRNAoutput.txt' or die "Can't open search file: $!";
@a=<$search_file2>;
for ($i=0;$i<=scalar @a;$i++){
if(@a[$i] =~ /^>/)
{
push(@possss,$i);

}
if(@a[$i] =~ m/^Query=/)
{
push(@Qpossss,$i);
}
}
$lastele=pop(@Qpossss);
push(@Qpossss,$lastele);
$lastele+=$lastele;
push(@Qpossss,$lastele);

for($i=0;$i<=(scalar @Qpossss);$i++){

for($j=0;$j<=scalar @possss;$j++)
{
if((@Qpossss[$i] < $possss[$j]) && ($possss[$j] < @Qpossss[$i+1]))
{
#print @Qpossss[$i]."\t".$possss[$j]."\n";
$x=@Qpossss[$i];
$a=@a[$x];
@B=split(" ",$a);
$B=">".@B[1];
$y=@possss[$j];
$length=@a[$y+4];
$length=substr($length,14,2);
if(($length > 17) && ($length < 26)){
$start=@a[$y+7];
$start=substr($start,7,2);
push (@keymatch,$B);
push (@keymatch,$start."\t".$length);
#print $B."\t".@possss[$j]."\t".$start."\t".$length."\n";
}
last;
}
}
}
%blast_miRNA=@keymatch;
open (READ,"< ".$path."\\miRNA\\pre-miRNA.fasta");
@data=<READ>;
%data=@data;
open(hh,"> ".$path."\\miRNA\\miRNA.txt");
foreach $key(keys %data)
{
$value=$data{$key};
$key=~s/^\s+|\s$//g;
foreach $blast_key(keys %blast_miRNA)
{
$blast_key2=$blast_key;
$blast_key2=~s/^\s+|\s$//g;
if ($blast_key2 eq $key){
$keymatch_value=$blast_miRNA{$blast_key};
#print $blast_key."\t".$keymatch_value."\n";

@position_align=split("\t",$keymatch_value);
$start=@position_align[0];
$end=@position_align[1];
$miRNA_pre=substr($value,$start-2,$end+2);
if(length($miRNA_pre)>19 && length($miRNA_pre)<26){
print hh $blast_key."\n";
#print $miRNA_pre."\t";
$miRNA_pre=~s/U/T/g;
print $miRNA_pre."\n";
print hh $miRNA_pre."\n";
}}}
}
}
#print @keymatch;
$cmd="perl post_miRNA.pl -i $folder.fasta";
system($cmd);
$cmd="perl filter_other_nc.pl -i $folder.fasta";
system($cmd);