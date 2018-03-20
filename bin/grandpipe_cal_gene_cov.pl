#!/usr/bin/perl -w
use strict;

die "$0 <BED>  <depth file>  <gene id file>  <RESULT>\n	depth file:\tfile from samtools depth\n" if $#ARGV!=3;
open BED2,"<$ARGV[0]" or warn "$!";
open OUT,">$ARGV[3]" or warn "$!";
my %hash_gene_sample;
my %hash_all_gene;
while(<BED2>){
	chomp;
	my $gene=(split /\./,(split /\s+/,$_)[-1])[0];
	$hash_all_gene{$gene}=1;
}
close BED2;
my%hash_geneid;
open GENEID,"<$ARGV[2]" or warn "$!";
while(<GENEID>){
	chomp;
	if($_=~/^#/){
		next;
	}
	my($geneid,$symbol,$synonyms)=(split)[1,2,4];
	my @gene_list=split /\|/,$synonyms;
	push @{$hash_geneid{$symbol}},$geneid;
	foreach my$ele(@gene_list){
		push @{$hash_geneid{$ele}},$geneid;
	}
}	
open BED,"<$ARGV[0]" or warn "$!";	
open FILE,"<$ARGV[1]" or warn "$!";
chomp(my $bed=<BED>);
chomp(my $file=<FILE>);
my %gene_length;
my %gene_cov;
my (%gene_cov5,%gene_cov10,%gene_cov20);
my $last_pos=0;
while(1){
		last if (!$bed or !$file);
		chomp($bed,$file);
		my($bed_chr,$bed_start,$bed_end,$bed_other_info)=split /\s+/,$bed;

		if($bed_chr eq "X"){
			$bed_chr=23;
		}elsif($bed_chr eq "Y"){
			$bed_chr=24;
		}else{
		}

		my $gene=(split /\./,$bed_other_info)[0];
		my($file_chr,$file_pos)=(split /\s+/,$file)[0,1];
		my $file_depth=(split /\s+/,$file)[-1];
#		print "$file_chr\t";

		if($file_chr eq "X"){
			$file_chr=23;
		}elsif($file_chr eq "Y"){
			$file_chr=24;
		}else{
		}

		$gene_length{$gene}+=$bed_end-$bed_start+1 if $last_pos ne $bed_start;
		$last_pos=$bed_start;
		if($bed_chr==$file_chr){
			if($file_pos>=$bed_start && $file_pos<=$bed_end){
				$gene_cov{$gene}++ if $file_depth>=1;
				$gene_cov5{$gene}++ if $file_depth>=5;
				$gene_cov10{$gene}++ if $file_depth>=10;
				$gene_cov20{$gene}++ if $file_depth>=20;
				$file=<FILE>;
			}elsif($file_pos<$bed_start){
				$file=<FILE>;
			}elsif($file_pos>$bed_end){
				$bed=<BED>;
			}else{
				die "error\n";
			}
		}elsif($bed_chr>$file_chr){
			$file=<FILE>;
		}elsif($bed_chr<$file_chr){
			$bed=<BED>;
		}else{
			die "error2\n";
		}
	
}
	foreach my$ele(keys %hash_all_gene){
		if(not exists $gene_length{$ele}){
			print "no coverage:\t$ele\n";
			$hash_gene_sample{$ele}=0;
			my $gene_id_list;
			if(not exists $hash_geneid{$ele}){
				$gene_id_list=-1;
			}else{
				foreach my$ele2(@{$hash_geneid{$ele}}){
					$gene_id_list.=$ele2.",";
				}
				$gene_id_list=~s/,$//;
			}
			print OUT "$gene_id_list\t$ele\t0\t0\t0\t0\n";
		}elsif(not exists $gene_cov{$ele}){
			print "No coverage:\t$ele\n";
			$hash_gene_sample{$ele}=0;
			my $gene_id_list;
			if(not exists $hash_geneid{$ele}){
				$gene_id_list=-1;
			}else{
				foreach my$ele2(@{$hash_geneid{$ele}}){
					$gene_id_list.=$ele2.",";
				}
				$gene_id_list=~s/,$//;
			}
			print OUT "$gene_id_list\t$ele\t0\t0\t0\t0\n";
		}else{
			my($len,$cov)=($gene_length{$ele},$gene_cov{$ele});
			my $cov5=$gene_cov5{$ele}?$gene_cov5{$ele}:0;
			my $cov10=$gene_cov10{$ele}?$gene_cov10{$ele}:0;
			my $cov20=$gene_cov20{$ele}?$gene_cov20{$ele}:0;
			my $tmp=$cov/$len;
			my $tmp5=$cov5/$len;
			my $tmp10=$cov10/$len;
			my $tmp20=$cov20/$len;
			my $gene_id_list;
			if(not exists $hash_geneid{$ele}){
				$gene_id_list=-1;
			}else{
				foreach my$ele2(@{$hash_geneid{$ele}}){
					$gene_id_list.=$ele2.",";
				}
				$gene_id_list=~s/,$//;
			}
			print OUT "$gene_id_list\t$ele\t$tmp\t$tmp5\t$tmp10\t$tmp20\n";
		}
	}
=pod
my $flag=0;
print "gene\\sample";
foreach my$gene(sort {$a cmp $b} keys %hash_gene_sample){
	print OUT "$gene";
	foreach my$sample(sort {$a cmp $b} keys %{$hash_gene_sample{$gene}}){
		print OUT "%.4f\t",$hash_gene_sample{$gene}{$sample};
		print "\t$sample" if $flag==0;
	}
	print OUT "\n";
	$flag=1;
}
print "\n";
=cut
