
### Converting 
Transfer counts files to a folder called counts and select the top hits for each sample.  

```
cd step_4_output

mkdir counts
cut -f2 F5_24.merged.RefSeq_annotated | sort | uniq -c > counts/control_F5T24.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_28.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F5T28.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_32.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F5T32.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_48.merged.RefSeq_annotated | sort | uniq -c > counts/control_F5T48.merged.RefSeq_annotated_BLAST.txt


cut -f2 F6_24.merged.RefSeq_annotated | sort | uniq -c > counts/control_F6T24.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_28.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F6T28.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_32.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F6T32.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_48.merged.RefSeq_annotated | sort | uniq -c > counts/control_F6T48.merged.RefSeq_annotated_BLAST.txt


cut -f2 F8_24.merged.RefSeq_annotated | sort | uniq -c > counts/control_F8T24.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_28.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F8T28.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_32.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F8T32.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_48.merged.RefSeq_annotated | sort | uniq -c > counts/control_F8T48.merged.RefSeq_annotated_BLAST.txt


###############
cut -f2 F5_52.merged.RefSeq_annotated | sort | uniq -c > counts/control_F5T52.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_176.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F5T176.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_240.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F5T240.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_248.merged.RefSeq_annotated | sort | uniq -c > counts/control_F5T248.merged.RefSeq_annotated_BLAST.txt


cut -f2 F6_52.merged.RefSeq_annotated | sort | uniq -c > counts/control_F6T52.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_176.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F6T176.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_240.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F6T240.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_248.merged.RefSeq_annotated | sort | uniq -c > counts/control_F6T248.merged.RefSeq_annotated_BLAST.txt


cut -f2 F8_52.merged.RefSeq_annotated | sort | uniq -c > counts/control_F8T52.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_176.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F8T176.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_240.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F8T240.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_248.merged.RefSeq_annotated | sort | uniq -c > counts/control_F8T248.merged.RefSeq_annotated_BLAST.txt

#################
cut -f2 F5_264.merged.RefSeq_annotated | sort | uniq -c > counts/control_F5T264.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_56.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F5T56.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_74_5.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F5T74A5.merged.RefSeq_annotated_BLAST.txt

cut -f2 F5_152.merged.RefSeq_annotated | sort | uniq -c > counts/control_F5T152.merged.RefSeq_annotated_BLAST.txt


cut -f2 F6_264.merged.RefSeq_annotated | sort | uniq -c > counts/control_F6T264.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_56.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F6T56.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_74_5.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F6T74A5.merged.RefSeq_annotated_BLAST.txt

cut -f2 F6_152.merged.RefSeq_annotated | sort | uniq -c > counts/control_F6T152.merged.RefSeq_annotated_BLAST.txt


cut -f2 F8_264.merged.RefSeq_annotated | sort | uniq -c > counts/control_F8T264.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_56.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F8T56.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_74_5.merged.RefSeq_annotated | sort | uniq -c > counts/experimental_F8T74A5.merged.RefSeq_annotated_BLAST.txt

cut -f2 F8_152.merged.RefSeq_annotated | sort | uniq -c > counts/control_F8T152.merged.RefSeq_annotated_BLAST.txt
```

The *merged.RefSeq_annotated_BLAST.txt file are the ones we analyzed in R.
These are provided in the parent __rna_seq__ folder of this repository
