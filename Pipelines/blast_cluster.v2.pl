#!/usr/bin/perl
use warnings;
use strict;

die "perl $0 [in.len.sort] [in.blast.cvg] [out.fasta] [%cvg] [%iden]\n" unless @ARGV == 5;
my ($len, $in_f, $out_f, $cov, $iden) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;
### print STDERR "Program $0 Start...\n";

my $similar = $iden;
my $coverage = $cov;

###
my (@seq, %len, %info) = ();
my $curr = 9999999;

open IN, $len or die $!;
while(<IN>){
	chomp;
	my @s = split /\s+/;
	
	die "$len MUST BE SORTED.\n" if $s[-1] > $curr;
	$curr = $s[-1];

	push @seq, $s[0];
	$len{$s[0]} = $s[1];
}
close IN;

###
open IN, $in_f or die $!;
while(<IN>){
	chomp;
	my @s = split /\s+/;
	
	next if $s[0] eq $s[1];
	next unless exists $len{$s[0]} and $len{$s[1]}; ###
	next unless $s[5] >= $similar;
	next unless $s[4] / $len{$s[1]} * 100 >= $coverage;
	next unless $len{$s[0]} >= $len{$s[1]};

	push @{$info{$s[0]}}, $_;
}
close IN;

###
my %has = %len;
my $num = 0;

open OT, ">$out_f" or die $!;
foreach my $i(0..$#seq){
	my $seq = $seq[$i];
	next unless exists $has{$seq};
	print OT ">Cluster $num\n"; $num++;
	print OT "0\t$len{$seq}nt, >$seq... *\n";

	next unless exists $info{$seq};
	my $a = 0;
	my @info = @{$info{$seq}};
	for(@info){
		my @s = split /\s+/;
		next unless exists $has{$s[1]};
		$a++;
		printf OT "$a\t$len{$s[1]}nt, >$s[1]... at %.2f/$s[5]\n", $s[4]/$len{$s[1]}*100;
		delete $has{$s[1]};
	}
}
close OT;


open IN, $out_f or die $!;
open OT, ">$out_f.list" or die $!; 
my $b = 0;
my $a = 0;
while(<IN>){
    if(!/^>/){
        chomp;
        my @s = split /\s+/;
        $_ =~ m/>(\S+)/ ;
        my $b = $& ;
        $b =~ s/\.\.\..*//;
        $b =~ s/>//;
        if($s[0] == 0){
            $a = $b;
        }
        printf OT "$a\t$b\n";
    }     
}
close IN;
close OT;

print STDERR "Program End...\n";
############################################################
sub function {
	return;
}
############################################################

