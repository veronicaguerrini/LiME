#!/usr/bin/perl -w

=head
 -----------------------------------------------------------------------------------------------------------------------
 The script needs only two files:
 - a FASTQ file whose name is saved in the variable input,
 - and a list of titles saved in the variable headings.
 The reads in the input fastq file are filtered in an ordered fashion, according to the list of titles stored in the headings file.
 
 REPLACE the -input- and -headings- file with your own file names.
 
 USAGE:
        perl filter_by_headers.pl > NEW_FASTQ.fq

 -----------------------------------------------------------------------------------------------------------------------
=cut

my ($input) = "setA2_1.fq" ;
my ($headings) = "title_setA2.txt" ;

my($buffer) = "";
open(FILE, $input ) or die("Error reading file, stopped");

my ($buffer1)="";
open(FILE1, $headings ) or die("Error list, stopped");
my ($len) =  `cat $headings | wc -l`;

#Read files
$buffer1 = <FILE1>;
$buffer = <FILE>;

my $counter =0;

while ($counter < $len){
    chomp $buffer1;
    if ("$buffer"=~ /\Q$buffer1/)    {
		#print "IS\n";
        print "$buffer";
        $buffer = readline( *FILE );
        print "$buffer";
        $buffer = readline( *FILE );
        print "$buffer";
        $buffer = readline( *FILE );
        print "$buffer";
		$buffer1=readline(*FILE1);
		$counter=$counter+1;
    }
    else{
        #print "IS NOT\n";
        $buffer = readline( *FILE );
        $buffer = readline( *FILE );
        $buffer = readline( *FILE );
    }
        $buffer=readline(*FILE);
}

close(FILE);
close(FILE1);





