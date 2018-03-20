#!/usr/bin/python
import sys
import re
if(len(sys.argv)!=4):
	if(len(sys.argv)>=2):
		if(sys.argv[1]=="-h" or sys.argv[1]=="-H"):
			print(sys.argv[0]+'\t<bed file>   <samtools depth file>   <output file>')
			exit(1)	
	print("Arguments error. Type \"-h/H\" to get the help")
	exit(1)
chr_length={"1":249250621,"2":243199373,"3":198022430,"4":191154276,"5":180915260,"6":171115067,"7":159138663,"8":146364022,"9":141213431,"10":135534747,"11":135006516,"12":133851895,"13":115169878,"14":107349540,"15":102531392,"16":90354753,"17":81195210,"18":78077248,"19":59128983,"20":63025520,"21":48129895,"22":51304566,"X":155270560,"Y":59373566}
if_file=open(sys.argv[1],'r')
#of_file=open(sys.argv[2],'w')
gap_array=[]
ss=1
ee=1
last_chr="0"
for line in if_file:
	line=line.strip()
	line_eles=line.split('\t')
	cur_chr_length=chr_length[line_eles[0]]
	start=int(line_eles[1])
	end=int(line_eles[2])
	if(last_chr==line_eles[0]):
		if(start>ss):
			ee=start-1
			gap_array.append(line_eles[0]+'\t'+str(ss)+'\t'+str(ee))
			ss=end+1
	else:
		if(last_chr!="0"):
			gap_array.append(last_chr+'\t'+str(ss)+'\t'+str(chr_length[last_chr]))
		ss=1
		ee=1
		last_chr=line_eles[0]
if_file.close()
gap_array.append(last_chr+'\t'+str(ss)+'\t'+str(chr_length[last_chr]))
#if_file=open(sys.argv[2],'r')
#of_file=open(sys.argv[2]+'2','w')
large_gap_array=[]
iter=0
ss=1
ee=1
last_chr="0"
last_start="0"
last_end="0"
for line in gap_array:
	line=line.strip()
	line_eles=line.split('\t')
	if(last_chr==line_eles[0]):
		iter+=1
		if(iter%10==1):
			ss=line_eles[1]
		if(iter%10==0):
			ee=line_eles[2]
			large_gap_array.append(line_eles[0]+'\t'+str(ss)+'\t'+str(ee))
	else:
		if(iter%10!=0):
			if(last_chr=="0"):
				ss=line_eles[1]
			else:
				large_gap_array.append(last_chr+'\t'+str(ss)+'\t'+last_end)
				iter=1
				ss=line_eles[1]
		last_chr=line_eles[0]
	last_start=line_eles[1]
	last_end=line_eles[2]
depth_file_if=open(sys.argv[2],'r')
index=0
line=depth_file_if.readline()
bed_file=open(sys.argv[1],'r')
bed_line=bed_file.readline()
of_file=open(sys.argv[3],'w')
iter=0
for index in range(0,len(large_gap_array)):
	iter+=1
	array_value=large_gap_array[index]
	if(bed_line.strip()=='' or line.strip()==''):
		break
	gap_eles=array_value.split('\t')
	gap_eles[0]=gap_eles[0].replace('X','23')
	gap_eles[0]=gap_eles[0].replace('Y','24')
	bindata={}
	binbed=[]
	while 1:
		if(bed_line.strip()==''):
			break;
		bed_line=bed_line.strip()
		bed_eles=bed_line.split('\t')
		bed_eles[0]=bed_eles[0].replace('X','23')
		bed_eles[0]=bed_eles[0].replace('Y','24')
		if(int(bed_eles[0])==int(gap_eles[0])):
			if(int(bed_eles[2])<int(gap_eles[2])):
				binbed.append(bed_line)
				bed_line=bed_file.readline()
			else:
				break
		elif(int(bed_eles[0])<int(gap_eles[0])):
			binbed=[]
			bed_line=bed_file.readline()
		else:
			break
	while 1:
		if(line.strip()==''):
			break;
		line=line.strip()
		line_eles=line.split('\t')
		line_eles[0]=line_eles[0].replace('X','23')
		line_eles[0]=line_eles[0].replace('Y','24')
		if(len(line_eles)<3):
			print(line)
		chromosome=line_eles[0]
		position=line_eles[1]
		depth=line_eles[2]
		if(int(chromosome)==int(gap_eles[0])):
			if(int(position)<int(gap_eles[2])):
				bindata[position]=depth
				line=depth_file_if.readline()
			else:
				break
		elif(int(chromosome)<int(gap_eles[0])):
			line=depth_file_if.readline()
		else:
			break
	for ele in binbed:
		binbed_eles=ele.split('\t')
		cov1,cov5,cov10,cov20,cov30=0,0,0,0,0
		for base in range(int(binbed_eles[1])+1,int(binbed_eles[2])+1):
			if(str(base) in bindata):
				depth=bindata[str(base)]
				if(int(depth)>=1):
					cov1+=1
				if(int(depth)>=5):
					cov5+=1
				if(int(depth)>=10):
					cov10+=1
				if(int(depth)>=20):
					cov20+=1
				if(int(depth)>=30):
					cov30+=1
		length=float(int(binbed_eles[2])-int(binbed_eles[1]))
		cov1_ratio=float(cov1)/length
		cov5_ratio=float(cov5)/length
		cov10_ratio=float(cov10)/length
		cov20_ratio=float(cov20)/length
		cov30_ratio=float(cov30)/length
		of_file.write(ele+'\t'+str(cov1_ratio)+'\t'+str(cov5_ratio)+'\t'+str(cov10_ratio)+'\t'+str(cov20_ratio)+'\t'+str(cov30_ratio)+'\n')
	bindata={}
if(bed_line.strip!=''):
	of_file.write(bed_line+'\t0\t0\t0\t0\t0\n')
	while(1):
		bed_line=bed_file.readline()
		bed_line=bed_line.strip()
		if(bed_line.strip()==''):
			break
		of_file.write(bed_line+'\t0\t0\t0\t0\t0\n')
else:
	for index2 in binbed:
		of_file.write(index2+'\t0\t0\t0\t0\t0\n')
of_file.close()