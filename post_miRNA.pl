#!/usr/bin/perl
use Getopt::Long;
my $inputfile;
my $parastr = &GetOptions("i=s{1}"=>\$inputfile);
my @folder=split(".fa",$inputfile);
my @folder2=split /\\/,$inputfile;
$folder_name=@folder2[0];
my $folder=@folder[0];
$raw_file=$folder."_other_ncRNA.fasta";
#print $raw_file;
$presor_file=$folder_name."\\miRNA\\_precusor.fasta";
#print $presor_file;
$miRNA_file=$folder_name."\\miRNA\\miRNA.txt";
#print $miRNA_file;
$output=$folder_name."\\miRNA\\miRNA_details_output_result.txt";
#print $output;
open(miRNA_file,"< $miRNA_file");
@miRNA_file=<miRNA_file>;
%miRNA=@miRNA_file;
open(raw,"< $raw_file");
@raw_data=<raw>;
%raw_dis=@raw_data;
open(he,"< $presor_file");
@pre_file_data=<he>;
$size=scalar @raw_data;
for($i=0;$i<$size;$i++)
{
push(@total_raw_seq_id,@raw_data[$i]);
#print @raw_data[$i]."\t".$i."\n";
$i++;
}
$raw_seq_count=scalar @total_raw_seq_id;
#print "\n total row seq : $raw_seq\n";
for($i=0;$i<$raw_seq_count;$i++)
{
$raw_seq=@total_raw_seq_id[$i];
$count=0;
$raw_seq=~s/^\s+|\s+$//g;
foreach $pre_file_data_line(@pre_file_data)
{
@pre_file_data_line=split("\t",$pre_file_data_line);
$pre_file_data_seqs=@pre_file_data_line[0];
@pre_file_data_seq=split("-",$pre_file_data_line);
$pre_file_data_seq=@pre_file_data_seq[0];
$pre_file_data_seq=~s/^\s+|\s+$//g;
if($raw_seq eq $pre_file_data_seq)
{
$count++;
}
}
if($count>0){
print $raw_seq."\t".$count."\n";
push(@seq_repeat,$raw_seq);
push(@seq_repeat,$count);
}
}
$count_pre_line=scalar @pre_file_data;
$repeat_seq_count=scalar @seq_repeat;
open(w,"> $output");
for($i=0;$i<$repeat_seq_count;$i++)
{ 
$repeat_seq_id=@seq_repeat[$i];
print w $repeat_seq_id."\t".@seq_repeat[$i+1]."\n";
foreach $key (%raw_dis)
{
$key2=$key;
$key2=~s/^\s+|\s+$//g;
$repeat_seq_id=~s/^\s+|\s+$//g;
if ($key2 eq $repeat_seq_id)
{
$value=$raw_dis{$key};
print w $value;
print w length($value)-1;
print w "\n";
}
}
print w "--"."\n";
for($j=0;$j<$count_pre_line;$j++)
{
$pre_file_data_line=@pre_file_data[$j];
@pre_file_data_line=split("\t",$pre_file_data_line);
$pre_file_data_seqs=@pre_file_data_line[0];
@pre_file_data_seq=split("-",$pre_file_data_line);
$pre_file_data_seq=@pre_file_data_seq[0];
$pre_file_data_seq =~ s/^\s+|\s+$//g;
if($repeat_seq_id eq $pre_file_data_seq)
{

print w @pre_file_data_line[1]."\t".@pre_file_data_line[2]."\t".@pre_file_data_line[3];
print w @pre_file_data[$j+1];
print w @pre_file_data[$j+2];
$y=0;
foreach $key(%miRNA)
{
$key2=$key;
$key2=~s/^\s+|\s+$//g;
$pre_file_data_seqs=~s/^\s+|\s+$//g;
if ($key2 eq $pre_file_data_seqs)
{
$value=$miRNA{$key};
$value=~s/^\s+|\s+$//g;
print w "$value";
$y++;
}
}
if ($y == 0){
print w "Not found";
}
print w "\n";
}
$j=$j+3;
}
$i++
}



