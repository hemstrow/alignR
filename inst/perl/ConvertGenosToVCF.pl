if ($#ARGV == 3) {
  $file = $ARGV[0];
  $outfile = $ARGV[1];
  $format = $ARGV[2];
  $bamlist = $ARGV[3];
}
else {
  die;
}

@vcf_meta = ();

open(FILE, "<$file") or die;
open(OUTFILE, ">$outfile") or die;
open(BAMFILE, "<$bamlist") or die;

# print header
$time=localtime;

print OUTFILE ("##fileformat=VCFv4.0\n" . "##fileDate=" . $time . "\n" . "##source=alignR_with_ANGSD\n" . "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
print OUTFILE ("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n" . "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
print OUTFILE ("#CHROM\t" . "POS\t" . "ID\t" "REF\t" . "ALT\t" . "QUAL\t" . "FILTER\t" . "INFO\t" . "FORMAT\t");

while(<BAMFILE>){
  chomp($_);
  print OUTFILE ("$_\t");
}
print OUTFILE "\n";

my $snp=1;
while(<FILE>){
  
  chomp($_);
  # print("$_");
  @loci = split(/\t/, $_);
  print "\n====================================\n";
  foreach $locus (@loci){
    print "$locus\t";
  }
  print("\n");
  
  # initialize and grab chromosome/position
  $ngeno = 0;
  @outline = ();
  $chrom = shift(@loci); 
  $pos = shift(@loci); 
  $maj = shift(@loci);
  $min = shift(@loci);
  push(@outline, $chrom); # #Chrom
  push(@outline, $pos); # POS
  push(@outline, "SNP_", .$snp); # ID
  $snp++;
  
  @outarray = ();
  $ac = 0;
  
  # for NN format, count genotypes and prepare output
  if($format eq "NN"){
    foreach $geno (@loci){
      if($geno ne "NN"){
        @alleles = split(//, $geno);
        
        $this_genotype_code = "";
        
        # allele 1
        if($alleles[0] eq $maj){
          $this_genotype_code = $this_genotype_code . "0/"
        }
        else{
          $ac++;
          $this_genotype_code = $this_genotype_code . "1/"
        }
        
        # allele 2
        if($alleles[1] eq $maj){
          $this_genotype_code = $this_genotype_code . "0"
        }
        else{
          $ac++;
          $this_genotype_code = $this_genotype_code . "1"
        }
        
        $ngeno++;
      }
      else{
        $this_genotype_code = "./.";
      }
      print("$this_genotype_code\t");
      
      push(@outarray, $this_genotype_code);
    }
  }
  
  # for numeric
  elsif($format eq "numeric"){
    foreach $geno (@loci){
      if($geno ne -1){
        if($geno eq 0){
          $this_genotype_code = "0/0";
        }
        elsif($geno eq 1){
          $this_genotype_code = "0/1";
          $ac++;
        }
        else{
          $this_genotype_code = "1/1";
          $ac += 2;
        }
        
        $ngeno++;
      }
      else{
        $this_genotype_code = "./.";
      }
      
      push(@outarray, $this_genotype_code);
    }
    
    $nal = $ngeno * 2;
    $prop_ac = $ac/$nal;
    # if 0 was actually the minor (which shouldn't happen, really...), flip
    $i = 0;
    if($prop_ac > 0.5){
      foreach $genocode (@outarray){
        $outarray[i] =~ s/1/0/ig;
        $outarray[i] =~ s/0/1/ig;
        $i++;
      }
    }
  }
  
  
  else{
    die("Unaccepted input format.\n");
  }
  
  # finalize array and print
  push(@outline, $maj); #REF
  push(@outline, $min); #ALT
  push(@outline, "."); #QUAL
  push(@outline, "PASS"); #FILTER
  $info_string = "NS=" . $ngeno . ";AC=" . $ac;
  push(@outline, $info_string);
  push(@outline, "GT");
  push(@outline, @outarray);
  
  print OUTFILE join("\t", @outline) . "\n";
}

close(OUTFILE); close(FILE);