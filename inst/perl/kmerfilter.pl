if ($#ARGV == 8){ # $# means that we take index of last element
  $library_file  = $ARGV[0];
  $fastq_fileRA = $ARGV[1];
  $fastq_fileRB = $ARGV[2];
  $k = $ARGV[3];#length of kmer
  $rare = $ARGV[4]; #rare kmer threshold #default = 1
  $abundant = $ARGV[5]; #abundant kmer threshold  #default = 20000
  $rare_count = $ARGV[6]; #rare in a row #default = .8 of kmer length
  $abundant_count = $ARGV[7]; #abundant count #default = .8 of kmers in a read
  $output = $ARGV[8];
  $paired = 1;
}
elsif ($#ARGV == 7){ # $# means that we take index of last element
  $library_file  = $ARGV[0];
  $fastq_fileRA = $ARGV[1];
  $k = $ARGV[2];#length of kmer
  $rare = $ARGV[3]; #rare kmer threshold #default = 1
  $abundant = $ARGV[4]; #abundant kmer threshold  #default = 20000
  $rare_count = $ARGV[5]; #rare in a row #default = .8 of kmer length
  $abundant_count = $ARGV[6]; #abundant count #default = .8 of kmers in a read
  $output = $ARGV[7];
  $paired = 0;
}
else {
  die("Not the right number of arguments. \n");
}


open(FILE, "<$library_file") or die("Could not locate library file. \n");
open(FILE1, "<$fastq_fileRA") or die("Could not locate fastq file. \n");
if($paired = 1){
  open(FILE2, "<$fastq_fileRB") or die("Could not locate fastq file. \n");
}
open(OUTFILE1, '>', $output ."_RA.fastq") or die("Could not open output file. \n");
if($paired = 1){
  open(OUTFILE2, '>', $output ."_RB.fastq") or die("Could not open output file. \n");
}

%kmerlibrary;

$throwaway = <FILE>;

while(<FILE>){
  chomp($_);
  @data = split(/ /, $_);
  $kmerlibrary{$data[0]} = $data[1];
}

#print($kmerlibrary{"AAATGCTTGACATGAATTTTCTT"});

while(<FILE1>){
  $header1RA = $_;
  $line2RA = <FILE1>;
  $line3RA = <FILE1>;
  $line4RA = <FILE1>;
  # print("========================= \n $line2RA \n");
  
  if($paired = 1){
    $header1RB = <FILE2>;
    $line2RB = <FILE2>;
    $line3RB = <FILE2>;
    $line4RB = <FILE2>;
  }
  
  $reject = 0;
  
  #check RA
  chomp($line2RA);
  @seqRA = split(//, $line2RA);
  $loop_len_RA = ($#seqRA) - ($k - 1);
  
  # print("Entering RA loop \n");
  $reject = check_fastq($loop_len_RA, @seqRA);
  
  if($reject == 1){
    # print("read rejected \n");
    next;
  }
  #else{print("RA accepted, moving to RB \n");}
  
  if($paired = 1){
    #check RB
    chomp($line2RB);
    @seqRB = split(//, $line2RB);
    $loop_len_RB = ($#seqRB) - ($k - 1);
  
    print("Entering RB loop \n");
    $reject = check_fastq($loop_len_RB, @seqRB);
    
    if($reject == 1){
      # print("read rejected \n");
      next;
    }
    # else{print("RB accepted \n");}
  }
  
  # print("printing reads to outfile \n");
  
  #if passes all filters in both RA and RB, write RA and RB reads to RA and RB output files
  print(OUTFILE1 $header1RA . $line2RA . "\n" .  $line3RA . $line4RA);
  
  if($paired = 1){
    print(OUTFILE2 $header1RB . $line2RB . "\n".  $line3RB . $line4RB);
  }
}

#subfunction: check fastq files
sub check_fastq{
  $loop_len = shift(@_);
  $rare_in_a_row = 0;
  $abundant_ = 0;
  $reject = 0;
  @seq = @_;
  for($i = 0;$i <= $loop_len ; $i++){
    
    #print("=========================================\n $i \n");
    
    $end_of_kmer = $i + $k - 1;
    $this_kmer = join("", @seq[$i..$end_of_kmer]);
    #print("$this_kmer \n");
    
    #query kmer library hash, get count of kmer
    $kmer_count = $kmerlibrary{$this_kmer};
    #print("$kmer_count \n");
    #if rare, add one to $rare_in_a_row
    if($kmer_count <= $rare){ 
      #print("This kmer is rare \n");      
      $rare_in_a_row++;
      #print("$rare_in_a_row \n");
    }
    else{
      $rare_in_a_row = 0;
    }
    #if abundant, add one to #$abundant_RA
    if($kmer_count >= $abundant){
     # print("This kmer is abundant \n");
      $abundant_++;
     # print("$abundant_ \n");
    }
    
    #if we trip filter (ie $rare_in_a_row_RA reaches 10), skip and move on
    if($rare_in_a_row >= $rare_count){
      $i = $loop_len + 1;
      $reject = 1;
      # print("read rejected. rare filter \n");
      next;
    }
    
    if($abundant_ >= $abundant_count){
      $i = $loop_len + 1;
      $reject = 1;
      # print("read rejected. abundant filter \n");
      next;
    }
    
    #print OUTFILE("$this_kmer\n"";)
  }
  return($reject);
}


close(OUTFILE1); close(OUTFILE2); close(FILE1); close(FILE2);