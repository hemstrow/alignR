if ($#ARGV == 3) {
  $file = $ARGV[0];
  $file2 = $ARGV[1];
  $output = $ARGV[2];
  $output2 = $ARGV[3];
} 
else {
  die;
}

open(FILE, "<$file") or die;
open(FILE2, "<$file2") or die;
open(OUTFILE, '>', $output) or die $!;
open(OUTFILE2, '>', $output2) or die $!;

# $i=1;

while (<FILE>) {
  # print("=================================================================\nread: $i\n\n");
  
  $header = $_;
  $line2 = <FILE>;
  $line3 = <FILE>;
  $line4 = <FILE>;
  
  $headerR2 = <FILE2>;
  $line2R2 = <FILE2>;
  $line3R2 = <FILE2>;
  $line4R2 = <FILE2>;
  

  @sheader1 = split(//, $header);
  @sheaderR2 = split(//, $headerR2);
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
    # print("Checking headers:\n");
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
    # print("\n\t equi: $equi\n");
    
    @sheader = split(/ /, $header);
    $tail = $sheader[$#sheader];
    
    # print("\n\t tail: $tail\n");
    
    # replace the common ending pattern 1:N:0:1 with /1 or /2
    if($tail =~ m/(\d):N:\d+:[\dA-Z]+\n/){
      $R = $1;
      # print("\t\tkey:$R\n");
      if($R eq 3 | $R eq 2){
        $header = $equi . "/" . 2 ."\n";
        $headerR2 = $equi . "/" . 1 ."\n";
      }
      elsif ($R eq 1){
        $header = $equi . "/" . 1 ."\n";
        $headerR2 = $equi . "/" . 2 ."\n";
      }
    }
    # print("New R1 header: $header\nNew R2 header: $headerR2;")
  }
  
  $i++;
  print OUTFILE $header . $line2 . $line3 . $line4;
  print OUTFILE2 $headerR2 . $line2R2 . $line3R2 . $line4R2;
  
}

close(FILE); close(OUTFILE); close(FILE2); close(OUTFILE2);