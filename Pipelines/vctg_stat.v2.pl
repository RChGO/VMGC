#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl $0 [*.bt.f] [*.taxinfo] [*.faa.ngenes] [*.bt.f.tax]\n" unless @ARGV == 4;
my ($in_f, $tax_f, $gene_f, $out_f) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;

my (%tax, %ngene, %info) = ();

open IN, $tax_f or die $!;
<IN>;
while(<IN>){ chomp; my @s = split /\t+/; $tax{$s[0]} = $s[-1];}
close IN;

open IN, $gene_f or die $!;
while(<IN>){ chomp; my @s = split /\t+/; $ngene{$s[0]} = $s[-1];}
close IN;

############################################################
open IN, $in_f or die $!;
while(<IN>){
	chomp;
	my @s = split /\t+/;
	next unless exists $tax{$s[1]};
	my $a = $tax{$s[1]};
	$info{"$s[0]\t$a"} = $s[2] unless exists $info{"$s[0]\t$a"};
}
close IN;

my (%s_num, %s_idt) = ();

open OT,">$out_f.gene" or die $!;
for(sort keys %info){
	my $id = $info{$_};
	print OT "$_\t$info{$_}\n"; ###
	s/^(\S+?)_\d+\t/$1\t/;
	$s_num{$_}++;
	$s_idt{$_}+=$id;
}
close OT;

open OT, ">$out_f" or die $!;
for(sort keys %s_num){
	my $a = $s_idt{$_}/$s_num{$_};
	my @s = split /\t+/;
	print OT "$_\t$s_num{$_} of $ngene{$s[0]}\t$a\n";
}
close OT;

print STDERR "Program End...\n";
############################################################

