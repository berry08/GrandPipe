#!/usr/bin/python
import sys
import re
depth_file_if=open(sys.argv[1],'r')
output=open(sys.argv[4],'w')
output.write('#gene\tgene_id\ttranscript\tcds_start\tcds_end\t1x\t5x\t10x\t20x\t30x\n')
geneid_file_if=open(sys.argv[3],'r')
gene_id={}
for line in geneid_file_if:
	line=line.strip()
	line_eles2=line.split('\t')
	gene_id[line_eles2[0]]=line_eles2[1]
geneid_file_if.close()
line=depth_file_if.readline()
for i in range(1,25):
	pos_depth={}
	while(1):
		if(line.strip()==''):
			print('read finished')
			break
		line=line.strip()
		line_eles=line.split('\t')
		line_chr=line_eles[0]
		if(line_eles[0]=='X'):
			line_chr='23'
		elif(line_eles[0]=='Y'):
			line_chr='24'
		else:
			pass
		if(line_chr==str(i)):
			pos_depth[line_eles[1]]=line_eles[2]
			line=depth_file_if.readline()
		else:
#			print(str(i)+'\t'+str(line_chr)+'\t'+line)
			break
	refgene_if=open(sys.argv[2],'r')
	for ref_line in refgene_if:
		ref_line=ref_line.strip()
		ref_line_eles=ref_line.split('\t')
		ref_line_eles[2]=ref_line_eles[2].replace('chr','')
		if(len(ref_line_eles[2])>2 or ref_line_eles[2]=='M'):
			continue
		if(ref_line_eles[2]=='X'):
			ref_line_eles[2]=23
		elif(ref_line_eles[2]=='Y'):
			ref_line_eles[2]=24
		else:
			ref_line_eles[2]=int(ref_line_eles[2])
		if(ref_line_eles[2]==i):
			total,cov1,cov5,cov10,cov20,cov30=0,0,0,0,0,0
			exon_starts=ref_line_eles[9].split(',')
			exon_ends=ref_line_eles[10].split(',')
			cds_start=ref_line_eles[6]
			cds_end=ref_line_eles[7]
			exon_number=ref_line_eles[8]
			cds_start_exon_index=-1
			cds_end_exon_index=-1
			for j in range(0,int(exon_number)):
				if(cds_start>=exon_starts[j] and cds_start<exon_ends[j]):
					cds_start_exon_index=j
					break
			for j in range(0,int(exon_number)):
				if(cds_end>exon_starts[j] and cds_end<=exon_ends[j]):
					cds_end_exon_index=j
					break
			if(cds_start_exon_index<0 or cds_end_exon_index<0):
				continue
			for j in range(cds_start_exon_index,cds_end_exon_index+1):
				cur_start=0
				cur_end=0
				if(j==cds_start_exon_index):
					cur_start=cds_start
				else:
					cur_start=exon_starts[j]
				if(j==cds_end_exon_index):
					cur_end=cds_end
				else:
					cur_end=exon_ends[j]
				total+=int(cur_end) - int(cur_start)
				if(int(cur_end)<int(cur_start)):
					break
				#print(str(ref_line_eles[2])+'\t'+cur_start+'\t'+cur_end+'\t'+ref_line_eles[12]+'\t'+ref_line_eles[1])
				for pos in range(int(cur_start)+1,int(cur_end)+1):
					if str(pos) in pos_depth:
						depth=int(pos_depth[str(pos)])
						if(depth>=1):
							cov1+=1
						if(depth>=5):
							cov5+=1
						if(depth>=10):
							cov10+=1
						if(depth>=20):
							cov20+=1
						if(depth>=30):
							cov30+=1
			if(total==0):
				continue			
			cov1_ratio=float(cov1)/float(total)
			cov5_ratio=float(cov5)/float(total)
			cov10_ratio=float(cov10)/float(total)
			cov20_ratio=float(cov20)/float(total)
			cov30_ratio=float(cov30)/float(total)
			geneid=''
			if(ref_line_eles[12] in gene_id):
				geneid=gene_id[ref_line_eles[12]]
			else:
				geneid='-1'
			output.write(ref_line_eles[12]+'\t'+geneid+'\t'+ref_line_eles[1]+'\t'+cds_start+'\t'+cds_end+'\t'+str(cov1_ratio)+'\t'+str(cov5_ratio)+'\t'+str(cov10_ratio)+'\t'+str(cov20_ratio)+'\t'+str(cov30_ratio)+'\n')
	refgene_if.close()
output.close()
