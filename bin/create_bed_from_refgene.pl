#!/usr/bin/perl -w
use strict;
open FILE,"<$ARGV[0]" or warn "$!";
open OUT,">$ARGV[1]" or warn "$!";
while(<FILE>){
	chomp;
	my($nmid,$chr,$cds_start,$cds_end,$exon_start,$exon_end,$gene)=(split)[1,2,6,7,9,10,12];
	if($cds_start eq $cds_end){
		next;
	}
	$chr=~s/^chr//;
	if($chr!~/^\d+$/ && ($chr ne "X") && ($chr ne "Y")){
		next;
	}else{
		$exon_start=~s/,$//;
		$exon_end=~s/,$//;
		my @starts=split /,/,$exon_start;
		my @ends=split /,/,$exon_end;
		my $flag=0;
		if((scalar @starts) ne (scalar @ends)){
			die "error\n";
		}else{
			if((scalar @starts) ==1){
				#die join('\t',@starts);
				$flag=1;
				print OUT "$chr\t$cds_start\t$cds_end\t$nmid\t$gene\n";
			}else{
				my $start_idx=(scalar @starts)-1;
				my $end_idx=(scalar @ends)-1;
				for(my $i=0;$i<scalar @starts;$i++){
					if($cds_start<$starts[$i]){
						$start_idx=$i-1;
						last;
					}
				}
				#die $start_idx,"\n" if $gene eq "PKD1";
				for(my $i=0;$i<scalar @ends;$i++){
					if($cds_end<=$ends[$i]){
						$end_idx=$i;
						last;
					}
				}
				#print "$start_idx\t$end_idx\n";
				for(my $i=$start_idx;$i<=$end_idx;$i++){
					print OUT "$chr\t";
					#print "$i\n";
					if($i==$start_idx){
						print OUT "$cds_start\t";
					}else{
						print OUT "$starts[$i]\t";
					}
					if($i==$end_idx){
						print OUT "$cds_end\t";
					}else{
						print OUT "$ends[$i]\t";
					}
					print OUT "$nmid\t$gene\n";
					$flag=1;
				}
			}
		}
		if($flag==0){
			die "here2\t$nmid\n";
		}
	}
}