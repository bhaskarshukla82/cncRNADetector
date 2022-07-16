
#$filename="other_nc_rna.txt";
$filename=<>;
chomp($filename);
open(r,"<","$filename");
@a=<r>;
%hash=@a;
foreach $key (keys %hash)
{
$value= $hash{$key};
@b=split("",$value);
$total= scalar @b;
$count=0;
foreach $b(@b)
{ 
  if(($b eq 'G')||($b eq 'g')||($b eq 'c')||($b eq 'C'))
{
   $count++; 
}}
  $gc= ($count/$total)*100;
   $gc{$key}=$gc; 
}
#print "  lncRNA ID   :\t GC%    \n";
#foreach $key (keys %gc)
#{
#$value=$gc{$key};
#chomp($key);
#print" $key  :  $value \n";
#}
open(h,">","$filename output.fasta");
#print h "  lncRNA ID   :\t GC%    \n";
foreach $key (keys %gc)
{
$gc=$gc{$key};
$seq=$hash{$key};
@bb=split("",$seq);
$totals= scalar @bb;
$totals=$totals-1;
$key=~ s/^\s+|\s+$//g;
print h "$key \t$gc \t$totals\n";

}

close(h);

