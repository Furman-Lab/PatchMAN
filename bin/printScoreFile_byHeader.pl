#!/usr/bin/perl

#index of $param1 in param2, with first index==0
sub indexArray{
 1 while $_[0] ne pop;
 $#_
}


#print the file $param1, and the columns labeled by @param2
sub printFile
  {
    my $fname=shift;
    open(INFILE, $fname) or die "couldn't open $fname\n";
    my $line;
    do {
	$line = <INFILE>;
    } while( ! ( $line =~ m/score/) );
    my @file_header=split(/\s+/, $line);
    my @indices;
    while($_=shift){
      push @indices, indexArray($_,@file_header);
    }
    my $line;
    while($line = <INFILE>){
      print "$fname ";
      chomp $line;  
      my @fields = split / +/,$line;
      for($i=0; $i <= $#indices ; $i++){
	print $fields[$indices[$i]]," ";
      }
      print "\n";
    }

  }

$file=shift;
# extract and print header
my @header;
print "src_file ";
while($_ = shift){
  push @header,$_;
  print "$_ ";
}
print "\n";
printFile($file,@header);

