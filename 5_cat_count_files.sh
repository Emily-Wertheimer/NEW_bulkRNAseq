#Add count files from each sample all together (on command line)

# find all files with '.txt' extension in wd; extract 3rd field; put extracted data into new file with '.raw' extention 
for i in `find . -name "*.txt" -type f`; do     awk {'print $3'} "$i" > "$i".raw; done       

# merge files horizontally: select all files of given naming convention and concatenate horizontally into new file
paste -d "\t" *count_matrix.txt.raw > RawCountMatrix.txt            

# use 1 sample to get gene symbols; put gene symbols in final raw count mat
awk {'print $1'} 5536C1B_168_025_S3.BAM_count_matrix.txt > GeneSymbolsMatrix.txt     

# creates a big matrix
paste -d "\t" *Matrix* > RawCountMatrix_final.txt            


# make column names not ridiculous (replace .bam with whatever you want to replace)
sed -i -e 's/.bam//g' RawCountMatrix_final.txt 
