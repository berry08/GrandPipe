cc=g++
all:GrandPipe cal_on_target_ratio filter_bam final_stat hdfilter merge_uni_hap_vcf process.indel.region stat.depth stat_homo
exe=GrandPipe
exe2=cal_on_target_ratio
exe3=filter_bam
exe4=final_stat
exe5=hdfilter
exe6=merge_uni_hap_vcf
exe7=process.indel.region
exe8=stat.depth
exe9=stat_homo
#exe10=simply_gender_decide
INCLUDE = .

obj=main.o gc.o PipeGraph.o create_pipe.o check_graph_exists.o pipeline.o global_cmd_para.o postprocess.o check.o run.o
$(exe):$(obj)
	$(cc) -o $(exe) $(obj)
main.o:main.cpp gc.h PipeGraph.h create_pipe.h check_graph_exists.h pipeline.h global_cmd_para.h postprocess.h
	$(cc) -c -O3 main.cpp

PipeGraph.o:PipeGraph.cpp check.h run.h gc.h
	$(cc) -c -O3 PipeGraph.cpp
create_pipe.o:create_pipe.cpp check.h run.h gc.h
	$(cc) -c -O3 create_pipe.cpp
check_graph_exists.o:check_graph_exists.cpp gc.h
	$(cc) -c -O3 check_graph_exists.cpp
pipeline.o:pipeline.cpp
	$(cc) -c -O3 pipeline.cpp
global_cmd_para.o:global_cmd_para.cpp gc.h
	$(cc) -c -O3 global_cmd_para.cpp
postprocess.o:postprocess.cpp global_cmd_para.h
	$(cc) -c -O3 postprocess.cpp
check.o:check.cpp global_cmd_para.h gc.h
	$(cc) -c -O3 check.cpp
run.o:run.cpp gc.h
	$(cc) -c -O3 run.cpp
gc.o:gc.cpp
	$(cc) -c -O3 gc.cpp

$(exe2):cal_on_target_ratio.o gc.o
	$(cc) -o $(exe2) cal_on_target_ratio.o gc.o
cal_on_target_ratio.o:cal_on_target_ratio.cpp gc.h
	$(cc) -c -O3 cal_on_target_ratio.cpp


$(exe3):filter_bam.o gc.o
	$(cc) -o $(exe3) filter_bam.o gc.o
filter_bam.o:filter_bam.cpp gc.h
	$(cc) -c -O3 filter_bam.cpp

$(exe4):final_stat.o gender_decide.o gc.o
	$(cc) -o $(exe4) final_stat.o gender_decide.o gc.o
gender_decide.o:gender_decide.cpp gc.h
	$(cc) -c -O3 gender_decide.cpp
final_stat.o:final_stat.cpp gender_decide.h gc.h
	$(cc) -c -O3 final_stat.cpp

$(exe5):hdfilter.o gc.o
	$(cc) -o $(exe5) hdfilter.o gc.o
hdfilter.o:hdfilter.cpp gc.h
	$(cc) -c -O3 hdfilter.cpp

$(exe6):merge_uni_hap_vcf.o gc.o
	$(cc) -o $(exe6) merge_uni_hap_vcf.o gc.o
merge.uni_hap.vcf.o:merge_uni_hap_vcf.cpp gc.h
	$(cc) -c -O3 merge_uni_hap_vcf.cpp

$(exe7): process.indel.region.o gc.o
	$(cc) -o $(exe7) process.indel.region.o gc.o
process.indel.region.o: process.indel.region.cpp gc.h
	$(cc) -c -O3 process.indel.region.cpp

$(exe8):stat.depth.o gc.o
	$(cc) -o $(exe8) stat.depth.o gc.o
stat.depth.o:stat.depth.cpp gc.h
	$(cc) -c -O3 stat.depth.cpp

$(exe9):stat.homo.sample.o gc.o
	$(cc) -o $(exe9) stat.homo.sample.o gc.o
stat.homo.sample.o:stat.homo.sample.cpp gc.h
	$(cc) -c -O3 stat.homo.sample.cpp
#$(exe10):

