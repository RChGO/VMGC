#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl $0 [*.hmm] [scores.cutoff] [*.hmm.out] [*.faa] [*.faa.busco]\n" unless @ARGV == 5;
my ($in_f, $score, $out_f, $faa, $fout) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;

my %score = ();
open IN, $score or die $!;
while(<IN>){
	chomp;
	my @s = split /\s+/;
	$score{$s[0]} = $s[1];
}
close IN;

my %busco = ();
open IN, $in_f or die $!;
open OT,">$out_f" or die $!;
while(<IN>){
	chomp;
	next if /^#/;
	my @s = split /\s+/;
	next unless exists $score{$s[2]} and $s[5] >= $score{$s[2]};
	print OT "$_\n";
	$busco{$s[0]}++;
}
close IN;
close OT;

my (%ctg, %ctgb) = ();
open IN, $faa or die $!;
while(<IN>){
	chomp;
	if(/^>(\S+?)\s/){
		my ($gene, $ctg) = ($1,$1);
		$ctg =~ s/_\d+$//;
		$ctg{$ctg}++;
		$ctgb{$ctg}++ if exists $busco{$gene};
	}
}
close IN;

open OT, ">$fout" or die $!;
foreach(keys %ctg){
	my $a = 0; $a=$ctgb{$_} if exists $ctgb{$_};
	print OT "$_\t$ctg{$_}\t$a\n";
}
close OT;

print STDERR "Program End...\n";
############################################################

