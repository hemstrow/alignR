if ($#ARGV == 1) {
  $file = $ARGV[0];
  $output = $ARGV[1];
} 
else {
  die;
}

open(FILE, "<$file") or die;
open(OUTFILE, '>', $output) or die $!;

while (<FILE>) {
  
  $header = $_;
  $line2 = <FILE>;
  $line3 = <FILE>;
  $line4 = <FILE>;
  
  @sheader = split(/ /, $header);
  $tail = $sheader[$#header];
  
  if($tail =~ m/(\d):N:\d:\d\n/){
    $R = $1;
    if($R eq 3){
      $R = 2;
    }
    pop(@sheader);
    push(@sheader, "/" . $R . "\n");
    
    $header = join("", @sheader);
  }
  
  print OUTFILE $header . $line2 . $line3 . $line4;
}

close(FILE); close(OUTFILE)