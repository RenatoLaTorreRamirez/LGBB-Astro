# LGBB-Astro
## Read Preprocessing
```

```
## Target capture analysis
### Prepare input using CDS as references
```
CDS_file="/path/to/fasta"

seqkit fx2tab $CDS_file | sed 's/-/\t/' | awk '{print $1 "-gene" "\t" (NR)}' |  awk '{printf($1 "%03d\n", $2 "\t" $4)}' > ${CDS_file%fasta}names
seqkit fx2tab $CDS_file | cut -f2 | paste ${CDS_file%fasta}names - | seqkit tab2fx > Names_$CDS_file
```
### Assemble reads into reference 'bubbles'
```
threads=30
species="/path/to/names"
reads="/path/to/reads_folder"

while read name; do
  hybpiper assemble -t_dna Named_$CDS_file -r $reads/${name}_trim-READ?.fastq --prefix $name --bwa --cpu $threads;
done < $species
```
### Gather assembly statistics
```
hybpiper stats -t_dna Named_$CDS_file gene $species --stats_filename hybpiper_stats_Heyduk --seq_lengths_filename seq_lengths_176
hybpiper recovery_heatmap seq_lengths_176.tsv --heatmap_filename recovery_heatmap_176
```
### Retrieve exon sequences
```
hybpiper retrieve_sequences dna -t_dna Named_$CDS_file --sample_names $species
mkdir exons
mv gene* exons/
```
## Target capture analysis using previous Hybpiper version (1.3)
### Prepare input using CDS as references
```
python Hybpiper1.3/HybPiper-1.3.1_final/reads_first.py --check-depend
```
### Map reads from each species to each reference gene
```
while read name; do
  python Hybpiper1.3/HybPiper-1.3.1_final/reads_first.py --bwa --cpu $threads -r $reads/${name}_trim-READ?.fastq -b Named_$CDS_file;
done < $species
```
### Retrieve intron sequences
```
python Hybpiper1.3/HybPiper-1.3.1_final/retrieve_sequences.py Named_$CDS_file $species dna
while read name; do
  python Hybpiper1.3/HybPiper-1.3.1_final/intronerate.py --prefix $name;
done < $species
python Hybpiper1.3/HybPiper-1.3.1_final/retrieve_sequences.py Named_$CDS_file . intron
python Hybpiper1.3/HybPiper-1.3.1_final/retrieve_sequences.py Named_$CDS_file --sample_names $samples intron
mkdir introns_old exons_old
mv gene*_introns.fasta introns_old/
mv gene* exons_old/
```
### Retrieve supercontig sequences
```
python Hybpiper1.3/HybPiper-1.3.1_final/retrieve_sequences.py Named_$CDS_file . supercontig
mkdir supercontigs_old
mv gene*_supercontig.fasta supercontigs_old/
```
## Trimming and tree analysis
```
# Put trimming commands
```
## BPP_analysis
```
mkdir Filtered_76ind_0Ns_10dels_1kb
cd Filtered_76ind_0Ns_10dels_1kb
sequences="/path/to/seq_folder"
table_names="/path/to/codes"

printf "Intron_Heyduk\nExon_Heyduk\nIntron_Louseau\nExon_Louseau" >> Datasets.txt

while read dataset; do
  mkdir -p $dataset;
  echo $dataset;
  for file in $(ls $sequences/$dataset/*.fasta | cut -d'/' -f11); do
    echo $file; while read line; do
      seq=$(echo $line | cut -d' ' -f2);
      name=$(echo $line | cut -d' ' -f1);
      grep -w $name $table_names | cut -f2 | sed "s/$/\t$seq/";
    done < <(cat $sequences/$dataset/$file | seqkit fx2tab) | seqkit tab2fx > $dataset/$file;
  done;
done < Datasets.txt

while read dataset; do
  echo $dataset;
  for file in $(ls $sequences/$dataset/*.fasta); do
    nseqs=$(cat $file | grep -c ">");
    length=$(cat $file | seqkit head -n 1 | seqkit stats --tabular | cut -f5 | sed -n '2p');
    dels=$(cat $file | seqkit head -n 1 | seqkit fx2tab | cut -f2 | grep -o "-" | wc -l)
    ns=$(cat $file | seqkit head -n 1 | seqkit fx2tab | cut -f2 | grep -o "n" | wc -l);
    printf "$file\t$nseqs\t$length\t$dels\t$ns\n";
  done > ${dataset}.info;
done < Datasets.txt

# Número de n = 0
# Gaps <= 0
# Número de seqs = 76
# Longitud de alineamiento >= 1000
for file in $(ls *.info);
  do
    cat $file | awk '$2 == 76' | awk '$5 == 0' | awk '$4 == 0' | awk '$3 >= 1000' | cut -f1;
  done > filter_files.txt

while read line; do
  source=$(echo $line | cut -d'/' -f3);
  ident=$(basename $line);
  cp $line ./${source}_$ident;
done < filter_files.txt

mkdir Phylip
cd Phylip
for i in $(ls ../*.fasta); do
  ./Fasta2Phylip.pl $i ${i%fasta}phylip;
  sed -i 's/^/\^/' ${i%fasta}phylip;
  tail -c +2 ${i%fasta}phylip >> BPP_in.txt;
  printf "\n" >> BPP_in.txt;
done
```
