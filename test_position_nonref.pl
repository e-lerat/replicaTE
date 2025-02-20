#!/usr/bin/perl
#use strict ;
use warnings ;
use Getopt::Long;

#program that compares the positions of non reference insertions detected by polymorphic TE detection programs and the expected positions of the non reference insertions in the deleted simulated chromosome produced by ReplicaTE

### initialisation des variables de fichiers
my ($in_file_annotation,$in_file_tool,$out_file1,$out_file2,$bp) ;

GetOptions('in1=s' => \$in_file_annotation, 'in2=s' => \$in_file_tool, 'out=s' => \$out_file, 'bp:i' => \$bp);


open FICHIER,"<$in_file_annotation" or die "Impossible to open $in_file_annotation.\n" ;
@ref=<FICHIER>;
close (FICHIER);

open FICHIER2,"<$in_file_tool" or die "Impossible to open $in_file_tool.\n" ;
@tool=<FICHIER2>;
close (FICHIER2);

open OUTFILE,">$out_file" or die "No $out_file\n" ;

#$n=-1;
$nb_total=$vrai_positif=$faux_positif=0;
$TE_previous="avant";
$start_previous = 0;
$dupli=$n = 0;

printf OUTFILE "TE\tcopy_nb\tTrue_Positif\tFalse_Positif\n";

foreach(@tool)
{
  @linetool=split(/,/, $_);
  chomp($linetool[$#linetool]);

  $TE_tool=$linetool[1];
  $n ++;
  if(($TE_previous ne $TE_tool) && ($TE_previous ne "avant"))
  {
    #print "$TE_tool\t$TE_previous\t$nb_total\t$vrai_positif\n";
    $faux_positif = $nb_total - $vrai_positif;
    printf OUTFILE "$TE_previous\t$nb_total\t$vrai_positif\t$faux_positif\n";
    $nb_total=$vrai_positif=$faux_positif=0;
  }
    $nb_total ++;
    #print "$nb_total\n";
    $start_tool=$linetool[3];
 
   if($TE_previous ne $TE_tool)
    {
     	foreach(@ref)
  	{
   		 @lineref=split(/[\t\_]/, $_);
    		chomp($lineref[$#lineref]);
    
   		 $TE=$lineref[0];
   		 $start=$lineref[2];
   		 
    		if($TE eq $TE_tool)
    		{
    			if(abs($start - $start_tool) < $bp) #start positions match with an error marging of x bp
    			{
    				$vrai_positif ++;
    			}
    		}
    	}
    	  $TE_previous=$TE_tool;
 	  $start_previous=$start_tool;	
   }
   
   elsif(($TE_previous eq $TE_tool) && (abs($start_tool - $start_previous >$bp)))
    {
     	foreach(@ref)
  	{
   		 @lineref=split(/[\t\_]/, $_);
    		chomp($lineref[$#lineref]);
    
   		 $TE=$lineref[0];
   		 $start=$lineref[2];
   		 
    		if($TE eq $TE_tool)
    		{
    			if(abs($start - $start_tool) < $bp) ##start positions match with an error marging of x bp
    			{
    				$vrai_positif ++;
    			}
    		}
    	}
    	  $TE_previous=$TE_tool;
 	  $start_previous=$start_tool;	
    
    }
    elsif(($TE_previous eq $TE_tool) && (abs($start_tool - $start_previous <$bp)))
    {
    	#print"$TE_previous\t$TE_tool\t$start_tool\t$start_previous\n";
    	$nb_total --;
    	$dupli ++;
    	#print "$nb_total\n";
    	$TE_previous=$TE_tool;
  	$start_previous=$start_tool;	

    }
 }

    $faux_positif = $nb_total - $vrai_positif;
    printf OUTFILE "$TE_previous\t$nb_total\t$vrai_positif\t$faux_positif\n"; 
    print "there was $dupli duplicated results in the output of the tools that were removed on a total of $n\n";

close(OUTFILE);
exit;
