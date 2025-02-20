#!/usr/bin/perl
#use strict ;
use warnings ;
use Getopt::Long;

#program that compares the positions of insertions detected by polymorphic TE detection programs and the positions of the reference insertions, making sure that the name of the copy is indeed in the complete simulated chromosome produced by ReplicaTE, to compute the number of True and False positives

### initialisation des variables de fichiers
my ($in_file_annotation,$in_file_tool,$in_file_annotation_del,$out_file1,$out_file2,$bp) ;

GetOptions('in1=s' => \$in_file_annotation, 'in2=s' => \$in_file_tool, 'in3=s' => \$in_file_annotation_del, 'out=s' => \$out_file, 'bp:i' => \$bp);


open FICHIER,"<$in_file_annotation" or die "Impossible to open $in_file_annotation.\n" ;
@ref=<FICHIER>;
close (FICHIER);

open FICHIER2,"<$in_file_tool" or die "Impossible to open $in_file_tool.\n" ;
@tool=<FICHIER2>;
close (FICHIER2);

open FICHIER3,"<$in_file_annotation_del" or die "Impossible to open $in_file_annotation_del.\n" ;
my @del=<FICHIER3>;
close (FICHIER3);
			
open OUTFILE,">$out_file" or die "No $out_file\n" ;

#$n=-1;
$nb_total=$vrai_positif=$faux_positif=0;
$TE_previous="avant";

printf OUTFILE "TE\tcopy_nb\tTrue_Positif\tFalse_Positif\n";

foreach(@tool)
{
  @linetool=split(/,/, $_);
  chomp($linetool[$#linetool]);

  $TE_tool=$linetool[1];
  
  if(($TE_previous ne $TE_tool) && ($TE_previous ne "avant"))
  {
    $faux_positif = $nb_total - $vrai_positif;
    printf OUTFILE "$TE_previous\t$nb_total\t$vrai_positif\t$faux_positif\n";
    $nb_total=$vrai_positif=$faux_positif=0;
  }
    $nb_total ++;
    
    $start_tool=$linetool[3];
    $end_tool=$linetool[4];
    
  foreach(@ref)
  {
    @lineref=split(/[\t\_]/, $_);
    chomp($lineref[$#lineref]);
    
    $TE=$lineref[0];
    $start=$lineref[2];
    $end=$lineref[3];


    if($TE eq $TE_tool)
    {
      #if((abs($start - $start_tool) < 20) && (abs($end - $end_tool) <20))
      if((abs($start - $start_tool) < $bp) && (abs($end - $end_tool) <$bp))

      {
      	foreach(@del)
      	{
      		my @linedel=split(/[\t]/, $_);
    		chomp($linedel[$#linedel]);
    		my $copy= $linedel[0];
    		my $name_ref=$TE . "_" . $lineref[1];

      		if($copy eq $name_ref)
      		{
      			$vrai_positif ++;
      		}	
      	}	
      }
    }
  }
  $TE_previous=$TE_tool;	

}
    $faux_positif = $nb_total - $vrai_positif;
    printf OUTFILE "$TE_previous\t$nb_total\t$vrai_positif\t$faux_positif\n"; 

close(OUTFILE);
exit;
