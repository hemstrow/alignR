if ($#ARGV == 3) {
  $fileR1 = $ARGV[0];
  $fileR2 = $ARGV[1];
  $barcode = $ARGV[2];
  $prefix = $ARGV[3];
}
else {
  die;
}

@commas = split(/\,/, $barcode);

# prepare hash for each barcode with the file names
$x=0;
while ($x <= $#commas){
  $hashR1{$commas[$x]} = $prefix . "_R1_" . $commas[$x] . ".fastq";
  $filenameR1 = $hashR1{$commas[$x]};
  open($filenameR1, ">$filenameR1") or die;
  
  $hashR2{$commas[$x]} = $prefix . "_R2_" . $commas[$x] . ".fastq";
  $filenameR2 = $hashR2{$commas[$x]};
  open($filenameR2, ">$filenameR2") or die;
  
  $x++;
}


open(FILER1, "<$fileR1") or die;
open(FILER2, "<$fileR2") or die;

# match barcodes to hash keys, print if a match
while (<FILER1>) {

  $R1_l1 = $_;
  $R1_l2 = <FILER1>;
  $R1_l3 = <FILER1>;
  $R1_l4 = <FILER1>;

  $R2_l1 = <FILER2>;
  $R2_l2 = <FILER2>;
  $R2_l3 = <FILER2>;
  $R2_l4 = <FILER2>;

  $bc = substr($R2_l2,0,6);

  if ($hashR1{$bc} ne ""){
    $F1 = $hashR1{$bc};
    $F2 = $hashR2{$bc};

    print $F1 $R1_l1 . $R1_l2 . $R1_l3 . $R1_l4;
    print $F2 $R2_l1 . $R2_l2 . $R2_l3 . $R2_l4;
  }
}


close FILER1; close FILER2;

$x=0;
while ($x <= $#commas){
  
  $F1 = $hashR1{$commas[$x]};
  $F2 = $hashR2{$commas[$x]};

  close($F1); close($F2);
  
  $x++;
}
