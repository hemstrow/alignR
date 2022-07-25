if ($#ARGV == 3) {
  $fileR1 = $ARGV[0];
  $fileR3 = $ARGV[1];
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
  
  $hashR3{$commas[$x]} = $prefix . "_R3_" . $commas[$x] . ".fastq";
  $filenameR3 = $hashR3{$commas[$x]};
  open($filenameR3, ">$filenameR3") or die;
  
  $x++;
}


open(FILER1, "<$fileR1") or die;
open(FILER3, "<$fileR3") or die;

# match barcodes to hash keys, print if a match
while (<FILER1>) {

  $R1_l1 = $_;
  $R1_l2 = <FILER1>;
  $R1_l3 = <FILER1>;
  $R1_l4 = <FILER1>;

  $R3_l1 = <FILER3>;
  $R3_l2 = <FILER3>;
  $R3_l3 = <FILER3>;
  $R3_l4 = <FILER3>;

  @bc = split(/:/, $R1_l1);
  $bc = @bc[$#bc]
  

  if ($hashR1{$bc} ne ""){
    $F1 = $hashR1{$bc};
    $F3 = $hashR3{$bc};

    print $F1 $R1_l1 . $R1_l2 . $R1_l3 . $R1_l4;
    print $F3 $R3_l1 . $R3_l2 . $R3_l3 . $R3_l4;
  }
}


close FILER1; close FILER3;

$x=0;
while ($x <= $#commas){

       $F1 = $hashR1{$commas[$x]};
       $F3 = $hashR3{$commas[$x]};

       close($F1); close($F3);

       $x++;
}
