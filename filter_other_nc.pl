#!/usr/bin/perl
use Getopt::Long;
my $inputfile;
my $parastr = &GetOptions("i=s{1}"=>\$inputfile);
my @folder=split(".fa",$inputfile);
my @folder2=split /\\/,$inputfile;
$folder_name=@folder2[0];
my $folder=@folder[0];
$raw_file=$folder_name."\\miRNA\\miRNA_details_output_result.txt";
$output=$folder_name."\\miRNA\\other_nc_after_miRNA.fasta";
$raw_file2=$folder."_other_ncRNA.fasta";
open(miRNA_file,"< $raw_file");
@miRNA_file=<miRNA_file>;
$miRNA_file=join("",@miRNA_file);
@miRNA_file2=split(">",$miRNA_file);
shift @miRNA_file2;
#print "@miRNA_file2[2]";
foreach $index(@miRNA_file2)
{
@id=split("\t",$index);
#print @id[0]."\n";
push(@miRNA_id,">".@id[0]);
}
open(other,"< $raw_file2");
@otherRNA_file=<other>;
%hash=@otherRNA_file;
$size=scalar @miRNA_id;
open(w,"> $output");
print $size;
foreach $key(keys %hash)
{
$key2=$key;
$key2=~ s/^\s+|\s+$//g;
if ( grep( /^$key2$/, @miRNA_id ) ) {

}
else{
$value=$hash{$key};
print w $key;
print w $value;
}
}

