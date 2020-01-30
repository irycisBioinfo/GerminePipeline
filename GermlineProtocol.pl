#!/usr/bin/perl

use Getopt::Long;
use strict;

my $r1;
my $r2;
my $ref;
my $aln;
my $vc;
my $platform;
my $prefix;
my $minLen;
my $minQual;
my $threads;
my $SeqRef;
my @Options;
my $dbSNP;
my $indexRef;
my $RG;
my $rmDups;
my $bed;
my $vcOpt;



@Options = (
		
		{OPT=>"r1=s",	VAR=>\$r1,	DESC=>"Reads 1"},
		{OPT=>"r2=s",	VAR=>\$r2,	DESC=>"Reads 2"},
		{OPT=>"ref=s",	VAR=>\$ref, DEFAULT => "GRCh37"	,DESC=>"Choose Reference genome (Ensembl -> GRCh37,NCBI -> GRCh38"},
		{OPT=>"aln=s", VAR=>\$aln, DEFAULT => "bwa", DESC=>"Aligner to use (bwa, bowtie2, novoaling)"},
		{OPT=>"vc=s", VAR=>\$vc, DEFAULT => "gatk", DESC=>"Choose the Variant Caller (gatk, varscan2, freebayes)"},
		{OPT=>"platform=s", VAR=>\$platform, DEFAULT=> "ILLUMINA", DESC=>"Specifie the sequencing platform (ILLUMINA,)"},
		{OPT=>"prefix=s", VAR=>\$prefix,DEFAULT =>"out", DESC=>"Prefix for output files"},
		{OPT=>"min_len=s", VAR=>\$minLen, DEFAULT => 140,DESC=>"Minimun read length"},
		{OPT=>"min_qual_mean=s", VAR=>\$minQual, DEFAULT => 20, DESC=>"Minimun read length"},
		{OPT=>"threads=s", VAR=>\$threads, DEFAULT => 20, DESC=>"Numbre of Threads to use"},
		{OPT=>"BED=s", VAR=>\$bed, DEFAULT => "", DESC=>"Use BED file to filter"},
		{OPT=>"duplicates=s", VAR=>\$rmDups, DEFAULT => "yes", DESC=>"Remove duplicates"},
		{OPT=>"vc_opt=s", VAR=>\$vcOpt, DEFAULT => "", DESC=>"Variant calling options (must be quoted '-option ...' "}
		
			);



#Check options and set variables
(@ARGV < 1) && (usage()); 
GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

# Now setup default values.file
foreach (@Options) {
	if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	${$_->{VAR}} = $_->{DEFAULT};
	}
}

if($ref eq "GRCh37")
{
	$SeqRef = "/storage/BD/eukaryotas/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa";
	$dbSNP = "/storage/BD/dbSNP/GRCh37/common_all_20180418.vcf";
	
	if($aln eq "bwa")
	{
		$indexRef ="/storage/BD/eukaryotas/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa";
		
	}elsif($aln eq "novoalign")
	{
		$indexRef = "/storage/BD/eukaryotas/Homo_sapiens/Ensembl/GRCh37/Sequence/NovoalingIndex/genome";
		
		
	}elsif($aln eq "bowtie2")
	{
		$indexRef = "/storage/BD/eukaryotas/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome";
		
	}
}elsif ($ref eq "GRCh38")
{
	$SeqRef = "/storage/BD/eukaryotas/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa";
    $dbSNP = "/storage/BD/dbSNP/GRCh38/common_all_20180418.vcf";
	if($aln eq "bwa")
	{
		$indexRef = "/storage/BD/eukaryotas/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa";
		
	}elsif($aln eq "novoalign")
	{
		$indexRef = "/storage/BD/eukaryotas/Homo_sapiens/NCBI/GRCh38/Sequence/NovoalingIndex/genome";
		
	}elsif($aln eq "bowtie2")
	{
		$indexRef = "/storage/BD/eukaryotas/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome";
		
	}
}else{
	print("\nReference genome not known\n");
	exit;
}



if($aln eq "bowtie2")
{
	$RG = " --rg-id $prefix --rg SM:$prefix --rg PL:$platform";
}else{
	$RG = "\@RG\\tID:$prefix\\tSM:$prefix\\tPL:$platform";
}

print("Filtering Reads ...\n");
#system("prinseq-lite.pl -fastq $r1 -fastq2 $r2 -min_len $minLen -min_qual_mean $minQual -out_good $prefix.good");
system("/storage/bioinfo/trimmomatic/trimmomatic-0.38.jar PE -threads $threads -phred33 $r1 $r2 $prefix.good_1.fastq bad_1 $prefix.good_2.fastq bad_2 SLIDINGWINDOW:15:$minQual MINLEN:$minLen");
system("rm *_prinseq_bad_*");

print("\n\nMapping reads...");
if($aln eq "bwa")
{
	system("bwa mem -R '$RG' -t $threads $indexRef $prefix.good_1.fastq $prefix.good_2.fastq > $prefix.$aln.sam");
}elsif($aln eq "novoalign")
{
	system("novoalign -f $prefix.good_1.fastq $prefix.good_2.fastq -o SAM '$RG' -d $indexRef -c $threads> $prefix.$aln.sam");
}elsif($aln eq "bowtie2")
{
    system("bowtie2 -x $indexRef -1 $prefix.good_1.fastq -2 $prefix.good_2.fastq -S $prefix.$aln.sam -p $threads $RG");
}

if( $bed eq "")
{
	system("samtools view --threads $threads -b -o $prefix.$aln.bam $prefix.$aln.sam ");
}else{
	system("samtools view -L $bed --threads $threads -b -o $prefix.$aln.bam $prefix.$aln.sam ");
}
system("samtools sort --threads $threads -o $prefix.$aln.sort.bam $prefix.$aln.bam");
if($rmDups eq "yes")
{
	system("gatk MarkDuplicates -I $prefix.$aln.sort.bam -M $prefix.metrix.dups -O $prefix.$aln.sort.dups.bam");
}else{
	system("cp $prefix.$aln.sort.bam $prefix.$aln.sort.dups.bam");
}
system("samtools index -\@ $threads $prefix.$aln.sort.dups.bam");
system("gatk BaseRecalibrator -I $prefix.$aln.sort.dups.bam -known-sites $dbSNP -output BQSR.table  -reference $SeqRef");
system("gatk ApplyBQSR --bqsr-recal-file BQSR.table -I $prefix.$aln.sort.dups.bam -O $prefix.$aln.sort.dups.bqsr.bam");

if($vc eq "gatk")
{
	system("gatk HaplotypeCaller -I $prefix.$aln.sort.dups.bqsr.bam -O $prefix.$aln.$vc.vcf -R $SeqRef $vcOpt");
}elsif( $vc eq "varscan2")
{
		system("samtools mpileup -B -f $SeqRef $prefix.$aln.sort.dups.bqsr.bam >$prefix.$aln.sort.dups.bqsr.pileup");
		system("java -jar /home/val/PerlBin/VarScan.v2.4.3.jar mpileup2cns $prefix.$aln.sort.dups.bqsr.pileup --variants --output-vcf 1 $vcOpt > $prefix.$aln.$vc.vcf");
		
}elsif($vc eq "freebayes")
{	
		system("freebayes -f $SeqRef $prefix.$aln.sort.dups.bqsr.bam $vcOpt > $prefix.$aln.$vc.vcf ");
}else{
	 print("Error: Unknown variant caller");
	 exit;
}


sub usage {
	foreach (@Options) {
		
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	print "\n\n\n";
	exit(1);
}
