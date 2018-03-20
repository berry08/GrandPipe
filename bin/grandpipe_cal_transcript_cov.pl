#!/usr/bin/perl -w
use strict;

die "$0 <BED>  <depth file>  <gene id file>  <RESULT>\n	depth file:\tfile from samtools depth\n" if $#ARGV!=3;
open BED2,"<$ARGV[0]" or warn "$!";
open OUT,">$ARGV[3]" or warn "$!";
my %hash_gene_sample;
my %hash_all_transcript;
while(<BED2>){
	chomp;
	my ($nmid,$gene)=(split /\s+/,$_)[3,4];
	$hash_all_transcript{$nmid}=$gene;
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
open FILE,"<$ARGV[1]" or warn "$!";
my @lines;
my %depth_hash;
my %transcript_length;
my (%transcript_cov,%transcript_cov5,%transcript_cov10,%transcript_cov20,%transcript_cov30);
my $iter=0;
my $length_flag=0;
while(<FILE>){
	chomp;
	my $info=$_;
	$iter++;
	if($iter%5000000!=0){
		my($chr,$pos,$depth)=split /\s+/,$info;
		$chr=$chr=~/X/?23:$chr=~/Y/?24:$chr;
		my $tmp_key=$chr*1000000000+$pos;
		$depth_hash{$tmp_key}=$depth;
	}else{
		#print "$_\n";
		open BED,"<$ARGV[0]" or warn "$!";	
		while(<BED>){
			chomp;
			my($bed_chr,$bed_start,$bed_end,$nmid,$gene)=split;
			$bed_chr=$bed_chr=~/X/?23:$bed_chr=~/Y/?24:$bed_chr;
			$transcript_length{$nmid}+=$bed_end-$bed_start if $length_flag==0;
			my $tmp_key2=$bed_chr*1000000000+$bed_start;
			my $tmp_key3=$bed_chr*1000000000+$bed_end;
			for(my $i=$tmp_key2+1;$i<=$tmp_key3;$i++){
				if(exists $depth_hash{$i}){
					if($depth_hash{$i}>=1){
						$transcript_cov{$nmid}++;
					}
					if($depth_hash{$i}>=5){
						$transcript_cov5{$nmid}++;
					}
					if($depth_hash{$i}>=10){
						$transcript_cov10{$nmid}++;
					}
					if($depth_hash{$i}>=20){
						$transcript_cov20{$nmid}++;
					}
					if($depth_hash{$i}>=30){
						$transcript_cov30{$nmid}++;
					}
				}
			}
		}
		$length_flag=1;
		close BED;
		undef %depth_hash;
		my($chr,$pos,$depth)=split /\s+/,$info;
		$chr=$chr=~/X/?23:$chr=~/Y/?24:$chr;
		my $tmp_key=$chr*1000000000+$pos;
		$depth_hash{$tmp_key}=$depth;
	}
}
open BED,"<$ARGV[0]" or warn "$!";	
while(<BED>){
	chomp;
	my($bed_chr,$bed_start,$bed_end,$nmid,$gene)=split;
	$bed_chr=$bed_chr=~/X/?23:$bed_chr=~/Y/?24:$bed_chr;

	$transcript_length{$nmid}+=$bed_end-$bed_start if $length_flag==0;
	my $tmp_key2=$bed_chr*1000000000+$bed_start;
	my $tmp_key3=$bed_chr*1000000000+$bed_end;
	for(my $i=$tmp_key2+1;$i<=$tmp_key3;$i++){
		if(exists $depth_hash{$i}){
			if($depth_hash{$i}>=1){
				$transcript_cov{$nmid}++;
			}
			if($depth_hash{$i}>=5){
				$transcript_cov5{$nmid}++;
			}
			if($depth_hash{$i}>=10){
				$transcript_cov10{$nmid}++;
			}
			if($depth_hash{$i}>=20){
				$transcript_cov20{$nmid}++;
			}
			if($depth_hash{$i}>=30){
				$transcript_cov30{$nmid}++;
			}
		}
	}
}
close BED;
$length_flag=1;
undef %depth_hash;


foreach my$ele(keys %hash_all_transcript){
	my $gene=$hash_all_transcript{$ele};
	if(not exists $transcript_length{$ele} || $transcript_length{$ele}==0){
		#print "no coverage:\t$ele\n";
		my $gene_id_list;
		if(not exists $hash_geneid{$gene}){
			$gene_id_list=-1;
		}else{
			foreach my$ele2(@{$hash_geneid{$gene}}){
				$gene_id_list.=$ele2.",";
			}
			$gene_id_list=~s/,$//;
		}
		print OUT "$gene_id_list\t$gene\t$ele\t0\t0\t0\t0\t0\n";
	}elsif(not exists $transcript_cov{$ele}){
		#print "No coverage:\t$ele\n";
		my $gene_id_list;
		if(not exists $hash_geneid{$gene}){
			$gene_id_list=-1;
		}else{
			foreach my$ele2(@{$hash_geneid{$gene}}){
				$gene_id_list.=$ele2.",";
			}
			$gene_id_list=~s/,$//;
		}
		print OUT "$gene_id_list\t$gene\t$ele\t0\t0\t0\t0\t0\n";
	}else{
		my($len,$cov)=($transcript_length{$ele},$transcript_cov{$ele});
		my $cov5=$transcript_cov5{$ele}?$transcript_cov5{$ele}:0;
		my $cov10=$transcript_cov10{$ele}?$transcript_cov10{$ele}:0;
		my $cov20=$transcript_cov20{$ele}?$transcript_cov20{$ele}:0;
		my $cov30=$transcript_cov30{$ele}?$transcript_cov30{$ele}:0;

		my $tmp=$cov/$len;
		my $tmp5=$cov5/$len;
		my $tmp10=$cov10/$len;
		my $tmp20=$cov20/$len;
		my $tmp30=$cov30/$len;
		$tmp=sprintf("%.4f",$tmp);
		$tmp5=sprintf("%.4f",$tmp5);
		$tmp10=sprintf("%.4f",$tmp10);
		$tmp20=sprintf("%.4f",$tmp20);
		$tmp30=sprintf("%.4f",$tmp30);
		if($tmp>1){
			$tmp=1;
		}
		if($tmp5>1){
			$tmp5=1;
		}
		if($tmp10>1){
			$tmp10=1;
		}
		if($tmp20>1){
			$tmp20=1;
		}
		if($tmp30>1){
			$tmp30=1;
		}
		my $gene_id_list;
		if(not exists $hash_geneid{$gene}){
			$gene_id_list=-1;
		}else{
			foreach my$ele2(@{$hash_geneid{$gene}}){
				$gene_id_list.=$ele2.",";
			}
			$gene_id_list=~s/,$//;
		}
		print OUT "$gene_id_list\t$gene\t$ele\t$tmp\t$tmp5\t$tmp10\t$tmp20\t$tmp30\n";
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
