#Get longest sequences for each UCE Probes
for i in $(ls -1 /path/to/uce/uce*nexus | cut -f1 -d ".")
do
Longest=$(sed 's/-//g' $i.nexus.fasta | esl-seqstat -a - | grep '=' | sort -k3,3n | tail -n1 | awk '{print $2}' )
Seq=$(grep $Longest $i.nexus | awk '{print $2}'| sed 's/\-//g')
printf ">$i\n" >>UCE_Probes.fasta
echo $Seq >>UCE_Probes.fasta
done

#Blast UCE onto genomes and get the position and number of hits blast version 2.13
for s in $(ls -1 /path/to/assemblies/*.fasta); 
do 
    Scaf=$(blastn -query UCE_Probes.fasta -db $s -num_threads 8 -outfmt 6 -evalue 1e-30 | awk '{if($10>$9){print $2,$9,$10,$10-$9} else{print $2,$10,$9,$9-$10}}' | sort -k4,4n | tail -n1 | awk '{print $1,$2,$3}'); 
    echo $i $Scaf >> $s.out;
		Scaf=$(blastn -query UCE_Probes.fasta -db $s -num_threads 8 -outfmt 6 -evalue 1e-30 | cut -f2 | sort -k1,1g | uniq | wc -l); 
	 	echo $i $Scaf >> $s.dup; 
done

for i in $(ls -1 *fasta.out)
do
	awk '!seen[$1]++' $i | awk 'BEGIN{OFS="\t"}{if($3-50<0){print $2,1,$4+50,$1}else{print $2,$3-50,$4+50,$1}}' > $i.bed
done

for s in $(ls -1 /path/to/assemblies/*.fasta | head -n7); 
do 
	blastn -query UCE_Probes.fasta -db $s -num_threads 2 -outfmt 6  -evalue 1e-30 | awk '{if($10>$9){print $1,$2,$9,$10,$10-$9} else{print $1,$2,$10,$9,$9-$10}}'| awk '{print $1,$2,$3,$4}' > $s.out 
done 

#Make a file with the count of each UCE per genome
echo UCE $(ls -1 *.fasta) >UCE.dup
awk '{print $1}' Biorhiza_pallida.fasta.dup > uce.dup
for i in $(ls -1 *fasta.dup)
do
	paste uce.dup $i | awk '{for(i=1;i<NF-1;i++){printf "%s ",$i};printf "%i\n",$NF}' > tmp.dup
	mv tmp.dup uce.dup 
done
cat uce.dup >> UCE.dup
rm uce.dup

#List UCE that have been duplicated in at least 4 genomes (potential paralogy)
awk '{if(NR>1){sum=0;for(j=2;j<=NF;j++){if($j>1)sum+=1};if(sum>3)print $1}}' UCE.dup > DuplicatedUCE.out

#generate bed and fasta file containing the UCEs -- need bedtools
mkdir UCE_seq_per_Species
for file in $(ls -1 /path/to/assemblies/*.fasta)
do
	cat DuplicatedUCE.out > UCEExclude.out
	awk '{if($2<1 || $2>1)print $1}' $file.dup >> UCEExclude.out
	grep -v -w -f UCEExclude.out $file.out | awk 'BEGIN{OFS="\t"}{if($3>0){print $2,$3,$4,$1}}' > $file\NoProbUCE.bed
	sed -i 's/\/path\/to\/uce_probes\/\///g' $file\NoProbUCE.bed
	sed -i 's/\.fasta//g' $file\NoProbUCE.bed
	bedtools getfasta -name -fi $file -bed $file\NoProbUCE.bed -fo UCE_seq/UCE$file
done

#Align UCEs with other species - need mafft
mkdir UCE_seq_per_UCE
for sp in $(ls -1 UCE_seq/UCE*fasta)
do
	Name=$(echo $sp | awk -F"/" '{print $NF}' | cut -f1 -d "." | sed 's/UCE//g')
	for uce in $(grep "^>" $sp  | sed 's/>//g')
	do
		printf ">$Name-Phylogenom\n" >> UCE_seq_per_UCE/$uce.nexus.fasta
		grep -A1 $uce $sp | tail -n1 >> UCE_seq_per_UCE/$uce.nexus.fasta
		done
done
mkdir Aligned_UCEs
for i in $(ls -1 UCE_seq_per_UCE/*.fasta); 
do 
	mafft --adjustdirection $i > Aligned_UCEs/$i; 
done
sed -i 's/_R_//g' Aligned_UCEs/uce*fasta

#Gblocks filtering of UCEs
for i in $(ls -1 Aligned_UCEs/uce-*fasta); 
do 
	./Gblocks $i -t n;
done

#Make nexus file for partition

for i in $(ls -1 Aligned_UCEs/*.htm); 
do 
	Name=$(echo $i | cut -f1 -d ".");
	NumPos=$(grep blue $i | cut -f2 -d "<" | sed 's/b>//g'); 
 	echo $Name $NumPos >> NumSitesGblocksFilt.txt; 
done

while read uce
do
	UCE=$(echo $uce | cut  -f1 -d ' ')
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $UCE.nexus.fasta-gb | sed 's/ //g' > Aligned_UCEs/$UCE.nexus.fasta
done <NumSitesGblocksFilt.txt

SpeciesList=$(grep "^>" Aligned_UCEs/*.fasta-gb | cut -f2 -d ":" | sort | uniq | sed 's/>//g')
for sp in $SpeciesList
do
echo ">"$sp > Concatenated_UCEs/$sp.fasta
while read uce
do
 UCE=$(echo $uce | cut  -f1 -d ' ')
 LENGTH=$(echo $uce | cut  -f2 -d ' ')
 MATCH=$(grep -c $sp Aligned_UCEs/$UCE.nexus.fasta-gb)
 if [[ $MATCH == 0 ]]
 then
   yes "n" | head -n $LENGTH | tr -d '\n' >> Concatenated_UCEs/$sp.fasta
 else
   awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Aligned_UCEs/$UCE.nexus.fasta-gb | grep -A1 $sp | tail -n1 | tr -d '\n' >> Concatenated_UCEs/$sp.fasta
 fi
done < NumSitesGblocksFilt.txt
done

NTAX=$(echo $SpeciesList | awk '{print NF}')
NCHAR=$(seqstat Acanthaegilips_brasiliensis_1339636.fasta | grep Largest | awk '{print $2}')
echo "BEGIN DATA;" > All_UCES.nexus
echo "DIMENSIONS  NTAX=$NTAX NCHAR=$NCHAR;" >> All_UCES.nexus
echo "FORMAT DATATYPE=DNA GAP=- MISSING=?;" >> All_UCES.nexus
echo "MATRIX" >> All_UCES.nexus
for sp in $SpeciesList
do
 Seq=$(tail -n1 $sp.fasta)
 printf "%s\t%s\n" $sp $Seq >> All_UCES.nexus
done
echo ";"  >> All_UCES.nexus
echo "" >> All_UCES.nexus
echo "END;" >> All_UCES.nexus
echo "" >> All_UCES.nexus
echo "BEGIN SETS;" >> All_UCES.nexus
echo "" >> All_UCES.nexus
printf "\t[loci]\n"  >> All_UCES.nexus
Start=0
End=0
count=1
while read uce
do
 UCE=$(echo $uce | cut  -f1 -d ' ')
 LENGTH=$(echo $uce | cut  -f2 -d ' ')
 Start=$(echo $End+1 | bc)
 End=$(echo $End+$LENGTH | bc)
 printf  "\tCHARSET %s = %s-%s;\n" $UCE $Start $End >> All_UCES.nexus
 charpart=$(echo $charpart $count:$UCE,)
 count=$(echo $count +1 | bc)
done < UCE_Length.txt
echo "" >> All_UCES.nexus
echo "	CHARPARTITION loci = $charpart" >> All_UCES.nexus
echo "END;" >> All_UCES.nexus

sed -i 's/,$/;/g' All_UCES.nexus
sed -i 's/uce-/uce_/g' All_UCES.nexus
sed -i 's/-Phy/_Phy/g' All_UCES.nexus

#Entropy-based automatic partitioning of UCE alignments -> https://github.com/Tagliacollo/PFinderUCE-SWSC-EN
python3 SWSCEN.py Conc/All_UCES.nexus




