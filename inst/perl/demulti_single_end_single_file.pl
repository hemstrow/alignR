if ($#ARGV == 4){
        $file1 = $ARGV[0];
        $barcode = $ARGV[2];
        $prefix = $ARGV[3];
        $barcode_in_header = $ARGV[4];
}
else {
        print "Incorrect number of arguments provided.\n";
        die;
}

@commas = split(/\,/, $barcode);
$barcode_length = length($commas[0]);

$x=0;
while ($x <= $#commas){
        $hash_r1{$commas[$x]} = $prefix . "_RA_" . $commas[$x] . ".fastq";
        $filename_r1 = $hash_r1{$commas[$x]};
        open($filename_r1, ">$filename_r1") or die;
        $x++;
}

open(FILE1, "<$file1") or die;


while (<FILE1>){

        $f1a = $_;
        $f1b = <FILE1>;
        $f1c = <FILE1>;
        $f1d = <FILE1>;
        
        if($barcode_in_header){
          @bc1 = split(/:/, $f1a);
          $bc1 = @bc1[$#bc1];
          
          out = $hash_r1{$bc1};
          print $out $f1a . $f1b . $f1c . $f1d;
        }
        else{
          $bc1 = substr($f1b,0,$barcode_length);
          $f1b_2 = substr($f1b, $barcode_length, length($f1b));
          $f1d_2 = substr($f1d, $barcode_length, length($f1d));
          
          $out = $hash_r1{$bc1};
          print $out $f1a . $f1b_2 . $f1c . $f1d_2;
        }
        
}
close FILE1;

$x=0;
while ($x <= $#commas){
        close($hash_r1{$commas[$x]});
        $x++;
}

