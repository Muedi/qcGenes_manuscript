from sys import *
import os
from multiprocessing import Process

assemblies = {'human': 'hg38', 'mouse': 'mm10', 'drosophila': 'dm6', 'celegans': 'ce10'}

bam_file_path = argv[1]
organism = argv[2]
output_file_ann = argv[3]
output_file_tss = argv[4]

# use the name of bam file for the temporary bed files
bam_file_basename = os.path.basename(bam_file_path)
mapping_bed_file_path = ('./%s.mapping.bed'%(bam_file_basename))
mapping_1M_bed_file_path = ('./%s.mapping_1M.bed'%(bam_file_basename))

#print('Converting bam to bed ...')
bedtools = ['bedtools', 'bamtobed', '-i', bam_file_path, '>', mapping_bed_file_path]
bedtools = ' '.join(bedtools)
os.system(bedtools)

#print('Sampling 1M reads from mapping.bed ...')
reads_sampling = ['shuf', '-n', '1000000', mapping_bed_file_path, '>', mapping_1M_bed_file_path]
reads_sampling = ' '.join(reads_sampling)
os.system(reads_sampling)

#run two Rscripts in parallel:

#print('\nRunning ChIPseeker')
annots_thread1 = ['Rscript', 'workflow/scripts/readsAnno.R', mapping_1M_bed_file_path, assemblies[organism], output_file_ann]
annots_thread1 = ' '.join(annots_thread1)
os.system(annots_thread1)

#print('\nRunning ChIPpeakAnno')
annots_thread2 = ['Rscript', 'workflow/scripts/readsAnnoTSS.R', mapping_1M_bed_file_path, assemblies[organism], output_file_tss]
annots_thread2 = ' '.join(annots_thread2)
os.system(annots_thread2)

#print('\nDeleting the temporary bed-files ...')
rm = ['rm', mapping_bed_file_path]
os.system(' '.join(rm))
rm = ['rm', mapping_1M_bed_file_path]
os.system(' '.join(rm))

