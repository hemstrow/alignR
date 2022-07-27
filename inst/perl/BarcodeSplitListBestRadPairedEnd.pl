if ($#ARGV == 4){
        $file1 = $ARGV[0];
        $file2 = $ARGV[1];
        $barcode = $ARGV[2];
        $prefix = $ARGV[3];
        $stacks_header = $ARGV[4];
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
        $hash_r2{$commas[$x]} = $prefix . "_RB_" . $commas[$x] . ".fastq";
        $filename_r1 = $hash_r1{$commas[$x]};
        $filename_r2 = $hash_r2{$commas[$x]};
        open($filename_r1, ">$filename_r1") or die;
        open($filename_r2, ">$filename_r2") or die;
        $x++;
}

open(FILE1, "<$file1") or die;
open(FILE2, "<$file2") or die;


while (<FILE1>){

        $f1a = $_;
        $f1b = <FILE1>;
        $f1c = <FILE1>;
        $f1d = <FILE1>;

        $f2a = <FILE2>;
        $f2b = <FILE2>;
        $f2c = <FILE2>;
        $f2d = <FILE2>;
        
        # fix headers to stacks style if requested
        if ($stacks_header){
          @sheader1 = split(//, $f1a);
          @sheaderR2 = split(//, $f2a);
          @equi_array = ();
  
          # check if we have to actually do anything
          $need_to_edit_headers = 1;
          if($sheader1[$#sheader1 - 1] eq "/" and $sheader1[$#sheader1] eq 1){
            if($sheaderR2[$#sheaderR2 - 1] eq "/" and $sheaderR2[$#sheaderR2] eq 2){
              $need_to_edit_headers = 0;
            }
          }
          elsif($sheader1[$#sheader1 - 1] eq "/" and $sheader1[$#sheader1] eq 2){
            if($sheaderR2[$#sheaderR2 - 1] eq "/" and $sheaderR2[$#sheaderR2] eq 1){
              $need_to_edit_headers = 0;
            }
          }
  
          if($need_to_edit_headers){
            # figure out matching header portion
            $prog = 0;
            foreach $h1p (@sheader1){
              if($h1p eq $sheaderR2[$prog]){
                if($h1p =~ /\S/){
                  push(@equi_array, $h1p);
                }
              }
              else{
                last;
              }
              $prog++;
            }

            $equi = join("", @equi_array);
    
            @sheader = split(/ /, $f1a);
            $tail = $sheader[$#sheader];
    
            # make the headers the common pattern + /1 or /2, easy since R1 is always /1
            $f1a = $equi . "/" . 1 ."\n";
            $f2a = $equi . "/" . 2 ."\n";
          }
        }

        # demultiplex
        $bc1 = substr($f1b,0,$barcode_length);
        $bc2 = substr($f2b,0,$barcode_length);

        if ($hash_r1{$bc1} ne "" && $hash_r1{$bc2} eq ""){

                $f1b_2 = substr($f1b, $barcode_length, length($f1b));
                $f1d_2 = substr($f1d, $barcode_length, length($f1d));

                $out1 = $hash_r1{$bc1};
                $out2 = $hash_r2{$bc1};

                print $out1 $f1a . $f1b_2 . $f1c . $f1d_2;
                print $out2 $f2a . $f2b . $f2c . $f2d;

        }
        elsif ($hash_r1{$bc1} eq "" && $hash_r1{$bc2} ne ""){

                $f2b_2 = substr($f2b, $barcode_length, length($f2b));
                $f2d_2 = substr($f2d, $barcode_length, length($f2d));

                $out1 = $hash_r1{$bc2};
                $out2 = $hash_r2{$bc2};

                print $out1 $f2a . $f2b_2 . $f2c . $f2d_2;
                print $out2 $f1a . $f1b . $f1c . $f1d;

        }
        elsif ($hash_r1{$bc1} ne "" && $hash_r1{$bc2} ne ""){

                print "Double Barcode!\t$bc1\t$bc2\n";

        }

}
close FILE1; close FILE2;

$x=0;
while ($x <= $#commas){
        close($hash_r1{$commas[$x]});
        close($hash_r2{$commas[$x]});
        $x++;
}

