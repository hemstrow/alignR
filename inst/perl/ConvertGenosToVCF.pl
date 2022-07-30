if ($#ARGV == 2) {
  $file = $ARGV[0];
  $outfile = $ARGV[1];
  $format = $ARGV[2];
}
else {
  die;
}

@vcf_meta = ();

open(FILE, "<$file") or die;
open(OUTFILE, ">$outfile") or die;

# print header
$time=localtime;

print OUTFILE ("##fileformat=VCFv4.0\n" . "##fileDate=" . $time . "\n" . "##source=alignR_with_ANGSD\n" . "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
print OUTFILE ("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele Count\">\n" . "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
print OUTFILE ("#CHROM\t" . "POS\t" . "REF\t" . "ALT\t" . "QUAL\t" . "FILTER\t" . "INFO\t" . "FORMAT\t")


while(<FILE>){
  
  $loci = chomp($_);
  @loci = split(/\t/, $loci);
  
  # initialize and grab chromosome/position
  $ngeno = 0;
  @outline = ();
  
  $chrom = shift(@loci);
  $pos = shift(@loci);
  push(@outline, $chrom);
  push(@outline, $pos);
  
  @outarray = ();
  
  # for NN format, count genotypes and prepare output
  if($format eq "NN"){
      %allele_hash = ();
      $a1 = "";
    
    foreach $geno (@loci){
      if($geno ne "NN"){
        @alleles = split(//, $geno);
        
        
        $this_genotype_code = "";
        foreach $allele (@alleles){
          if($a1 eq ""){
            $a1 = $allele;
          }
          
          if($a1 eq $allele){
            $this_genotype_code = $this_genotype_code . "A";
          }
          else{
            $this_genotype_code = $this_genotype_code . "B";
          }
          $allele_hash{$allele}++;
        }     
        
        
        $ngeno++;
        
        @this_genotype_code = split(//, $this_genotype_code);
        $this_genotype_code = join("/", @this_genotype_code);
      }
      else{
        $this_genotype_code = "./.";
      }
      
      push(@outarray, $this_genotype_code);
    }
    
    # figure out major and minor alleles
    @ordered_alleles = sort { $allele_hash{$a} <=> $allele_hash{$b} } keys %allele_hash;
    $maj = $ordered_alleles[2]; # major, is REF
    $min = $ordered_alleles[1]; # minor, is ALT
    
    # figure out ac
    $ac = $allele_hash{$min};
    
    # fix outarray
    $i <- 0;
    if($a1 eq $maj){
      foreach $genocode (@outarray){
        $outarray[i] =~ s/A/0/ig;
        $outarray[i] =~ s/B/1/ig;
        $i++;
      }
    }
    else{
      foreach $genocode (@outarray){
        $outarray[i] =~ s/A/1/ig;
        $outarray[i] =~ s/B/0/ig;
        $i++;
      }
    }
  }
  
  # for numeric
  elsif($format eq "numeric"){
    
    $ac = 0;
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
    
    
    $maj="0";
    $min="1";
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
  push(@outline, "GT")
  push(@outline, $outarray);
  
  print OUTFILE join("\t", @outline)
}

close(OUTFILE); close(FILE);