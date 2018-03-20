#ifndef CALLVCFFROMBAM_H
#define CALLVCFFROMBAM_H

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include "gc.h"
#include <unistd.h>
#include <stdlib.h>
using namespace::std;
void samtools_index(string root_dir,string bam_file);
void callVCF(string root_dir,string sample_name,string bam_file,string vcf_file,...);

#endif