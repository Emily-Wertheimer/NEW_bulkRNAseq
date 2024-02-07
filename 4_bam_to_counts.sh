#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --mem=200g
#SBATCH --j bam_to_counts # job name
#SBATCH --array=1-77 # number of participants as range starting at 1 (i.e., for 5 participants: 1-5)
#SBATCH --time=23:00:00 # HPC will give you this amount of time to run the process. This is usually enough time
#SBATCH -n 1 # how many nodes you are asking. This is running each subject on a differnt node so 1 is enough
#SBATCH --cpus-per-task=10 # How many CPUs. This is enough cpus no need to change

# resouce you are using are nodes * CPUs * memory - if you go above 120 per sample you will have to wait a lot of time to get an opening

# Outputs ----------------------------------
#SBATCH -o %x-%A-%a.out # this will give you the list of commands and there results (success/failure). If the run fails here you will get the spesifcs
#SBATCH -e %x-%A-%a.err # this will give you a short file with what errors were during the execution
#SBATCH --mail-user=emily.wertheimer@yale.edu # replace with your email
#SBATCH --mail-type=ALL

# ------------------------------------------
# put sample IDs here: 
sample_id=(31821A1_106_087_S24.BAM
31821B3_164_029_S35.BAM 
31821C1_130_063_S22.BAM 31822A2_141_052_S29.BAM 31822B1_129_064_S30.BAM 31822B4_154_039_S20.BAM 31823A1_117_076_S31.BAM 31823A3_189.BAM 31824A2_177_016_S26.BAM 31824A4_176_017_S34.BAM 
31824C2_118_075_S23.BAM  31824C3_188_005_S33.BAM 
31825B1_153_040_S28.BAM  31825B3_165_028_S27.BAM 
31825C3_142_051_S21.BAM  31826B2_152_041_S36.BAM 
31826C2_166_027_S19.BAM  31826C4_105_088_S32.BAM 
5531A1A_119_074_S15.BAM  5531C3C_190_003_S17.BAM 
5531C4B_132.BAM          5532B2B_107_086_S16.BAM 
5532B3C_143_050_S13.BAM  5532B4A_120_073_S7.BAM  
5533A2B_325_252_S12.BAM 5533A3A_167_026_S11.BAM 
5533C1A_179_014_S10.BAM 5534A1B_108_085_S8.BAM  
5534A3B_131_062_S14.BAM  5534C2C_178_015_S18.BAM 
5535B1A_337_240_S9.BAM 5535B3B_156_037_S4.BAM  
5535B4C_349_228_S2.BAM 5536C1B_168_025_S3.BAM  
5536C3B_144_049_S5.BAM 5536C4A_361_216_S1.BAM  
5As1A1A_272_305_S38.BAM 5As1A2B_211_366_S51.BAM 
5As1B3A_212_365_S43.BAM 5As2A3B_199_378_S52.BAM 
5As2B1A_224_353_S42.BAM 5As2C3A_247_330_S48.BAM 
5As3A4A_236_341_S41.BAM 5As3A4B_282_295_S53.BAM 
5As3C1A_283_294_S45.BAM 5As4A2A_260_317_S39.BAM 
5As4B1B_270_307_S54.BAM 5As4B4A_200_377_S44.BAM 
5As5B2B_284_293_S37.BAM 5As5C2A_259_318_S47.BAM 
5As5C4A_235_342_S49.BAM 5As6A1B_223_354_S50.BAM 
5As6A3A_248_329_S40.BAM 5As6B3B_271_306_S46.BAM 
NGN2A1A1_128_065_S38.BAM NGN2A1B3_175_018_S42.BAM
NGN2A2A2_174_019_S50.BAM NGN2A2B4_173_020_S58.BAM
NGN2A3A3_114_079_S55.BAM NGN2A3C1_150_043_S52.BAM
NGN2A4A4_115_078_S47.BAM NGN2A4C2_102_091_S56.BAM
NGN2A5B1_162_031_S51.BAM NGN2A5C3_313_264_S57.BAM
NGN2A6B2_116_077_S39.BAM NGN2A6C4_186_007_S49.BAM
NGN2B1A3_161_032_S59.BAM NGN2B1C1_127_066_S46.BAM
NGN2B2B1_138_055_S53.BAM NGN2B2B3_104_089_S40.BAM
NGN2B3A1_163_030_S43.BAM NGN2B3C3_140_053_S37.BAM
NGN2B4A2_151_042_S44.BAM NGN2B5B2_126_067_S54.BAM
NGN2B5B4_139_054_S45.BAM NGN2B6C2_187_006_S41.BAM
NGN2B6C4_103_090_S48.BAM)

# Directory containing BAM files
BAM_DIR="/gpfs/gibbs/pi/huckins/ekw28/running_alignment_jan27/bam_files/"

# this is where the magic starts
echo Starting ${sample_id[$SLURM_ARRAY_TASK_ID-1]}

# this is the line that runs the code. 
module load R/4.3.0-foss-2020b 
Rscript /gpfs/gibbs/pi/huckins/ekw28/counts/scripts/bam_to_counts.R ${sample_id[$SLURM_ARRAY_TASK_ID-1]}
