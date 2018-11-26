#!/usr/bin/python3
import sys
import os
software = "SAVE-ME v.0.09"
print("Welcome to {} !\nI'll now begin checking various files you've fed me.\nAnd will perform any additional steps, if required.".format(software))
citations="[1]. {}:\nWriting....\n[2]. GATK:\nThe Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data McKenna A,\nHanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA,\n2010 GENOME RESEARCH 20:1297-303\n".format(software)
files=[]
###############################################################################         EDIT THESE LINES            ########################################################################################################################################################################################################################################################################
#----------------------------------------Files------------------------------------------------#
threads=0       # '0' Will try to utilize all available CPU cores where it can. Make sure OpenMP is implemented properly on your system.
gatk = 'java -Xmx120G -jar /data/aditya/softwares/gatk-4.0.9.0/gatk-package-4.0.9.0-local.jar'
samtools = "samtools"
genome_assembly = "/data/aditya/softwares/gatk-4.0.9.0/datasets/GRCh38.fasta" #This SHOULD end in FASTA!
dbsnp = "/data/aditya/softwares/gatk-4.0.9.0/datasets/dbsnp_146.hg38.vcf"
bed_file= "/data/aditya/softwares/IBMFS.bed"      #Target regions to scan. Greatly speeds up the process and avoids unwanted alignments. Required. You can also use a Text file instead. Read: https://software.broadinstitute.org/gatk/documentation/article?id=11009
filter_snp= '--filter-expression "QD < 2.0" --filter-name "QD_FAIL" --filter-expression "FS > 60.0" --filter-name "FS_FAIL" --filter-expression "MQ < 40.0" --filter-name "MQ_FAIL" --filter-expression "ReadPosRankSum < -1.40" --filter-name "RPRS_FAIL" --filter-expression "SOR >=3.0" --filter-name "SOR_FAIL" --filter-expression "BaseQRankSum < -1.40" --filter-name "BQRS_FAIL" --filter-expression "DP <= 30.0" --filter-name "DP_FAIL" --filter-expression "MQRankSum <= -0.05" --filter-name "MQRS_FAIL"'
filter_indel='--filter-expression "QD < 2.0" --filter-name "QD_FAIL" --filter-expression "FS > 60.0" --filter-name "FS_FAIL" --filter-expression "MQ < 40.0" --filter-name "MQ_FAIL" --filter-expression "ReadPosRankSum < -1.55" --filter-name "RPRS_FAIL" --filter-expression "SOR >=3.0" --filter-name "SOR_FAIL" --filter-expression "BaseQRankSum < -1.00" --filter-name "BQRS_FAIL" --filter-expression "DP <= 30.0" --filter-name "DP_FAIL" --filter-expression "MQRankSum <= -0.05" --filter-name "MQRS_FAIL"'
effector ="java -Xmx120G -jar /data/aditya/softwares/snpeff/clinEff/ClinEff.jar -v -c /data/aditya/softwares/snpeff/clinEff/clinEff.config GRCh38.86 "
#effector ="java -Xmx20G -jar /Users/aditya/softwares/snpEff/snpEff.jar -v -c /Users/aditya/softwares/snpEff/snpEff.config GRCh38.86 "
readgroup="'@RG\\tID:MedGenome\\tCN:IlluminaHiSeq\\tSM:MedExome\\tPG:{}'".format(software)
bwa="bwa"
picard = "picard-tools" #Suggestion: Increase the amount of available RAM to picard. If PICARD is in your $PATH, edit it in any of your favorite text editor and add "-Xmx20G" for 20GB RAM or as you seem suitable.
gunzip="gunzip"
md5sum="md5sum"
fastqc="fastqc "
###############################################################################         DO NOT EDIT BELOW THIS      #########################################################################################################################################################################################################################################################################

#----------------------------------------Check Assembly Index--------------------------------------------------------------------------------------------------------------------------------#
print("*******************************\t\tInternal Tests Initiated\t************************************************************")
if os.path.isfile(genome_assembly):
    print("The Reference Genome Assembly '{}' seems available.\nProceeding.....\n".format(genome_assembly))
    pass
else:
    print("The Reference Genome Assembly '{}' isn't there where you pointed me!\nExiting...".format(genome_assembly))
    exit()
if os.path.isfile(genome_assembly+".fai"):
    print("The Reference Genome Assembly '{}' seems indexed.\nProceeding.....\n".format(genome_assembly))
    pass
else:
    print("The Reference Genome Assembly '{}' doen't seem to be indexed.\nSummoning SAMTOOLS!\nProceeding....".format(genome_assembly))
    os.system(samtools+" faidx "+genome_assembly)
    if os.path.isfile(genome_assembly + ".fai"):
        print("FASTA index file generated. This had to happen only once, provided you don't go and delete that INDEX file.")
    else:
        print("Well, this is embarassing!\nSome crazy thing happened.\nIndex file generation failed.\nPlease check if you have installed and correctly configured SAMTOOLS.\nAborting.")
        exit()
if os.path.isfile(genome_assembly+".ann") and os.path.isfile(genome_assembly+".sa") and os.path.isfile(genome_assembly+".amb") and os.path.isfile(genome_assembly+".pac") and os.path.isfile(genome_assembly+".bwt"):
    print("The Reference Genome Assembly '{}' seems BWA indexed.\nProceeding.....\n".format(genome_assembly))
    pass
else:
    print("The Reference Genome Assembly '{}' doen't seem to be indexed by BWA.\nSummoning BWA!\nProceeding....".format(genome_assembly))
    os.system(bwa+" index "+genome_assembly)
    if os.path.isfile(genome_assembly+".ann") and os.path.isfile(genome_assembly+".sa") and os.path.isfile(genome_assembly+".amb") and os.path.isfile(genome_assembly+".pac") and os.path.isfile(genome_assembly+".bwt"):
        print("FASTA BWA index file generated. This had to happen only once, provided you don't go and delete that INDEX file.\nI'll not recommend that.")
    else:
        print("Well, this is embarassing!\nSome crazy thing happened.\nBWA Index file generation failed.\nPlease check if you have installed and correctly configured BWA.\nAborting.")
        exit()

if os.path.isfile(genome_assembly[:-5]+"dict"):
    print("The Reference Genome Assembly '{}' seems to have a valid dictionary.\nProceeding.....\n".format(genome_assembly))
    pass
else:
    print("The Reference Genome Assembly '{}' doen't seem to have a valid dictionary.\nSummoning PICARD-TOOLS!\nProceeding....".format(genome_assembly))
    os.system(picard+" R="+genome_assembly+" O="+genome_assembly[:-5]+"dict")
    if os.path.isfile(genome_assembly[:-5] + "dict"):
        print("Dictionary file created.\nThis had to happen only once, provided you don't go and delete that dictionary file.")
    else:
        print("Well, this is embarassing!\nSome crazy thing happened.\nDictioinary file generation failed.\nPlease check if you have installed and correctly configured PICARD-TOOLS.\nAborting.")
        exit()
print("*******************************\t\tInternal Tests Completed\t************************************************************")


#----------------------------------------Assign Threads-----------------------------------------------------------------------------------------------------------------------------------------#

if threads == 0:
    threads = os.cpu_count()
else:
    pass

#==============================================MODULES==========================================================================================================================================#
                    # Thou shal PASS
def only_pass(vcf):
    print("Cleaning the variants which failed our filters.\nFile:{}".format(vcf))
    outfile=open(vcf[:-4]+"_PASS.vcf",'a')
    with open(vcf, 'r') as file:
        for line in file.readlines():
            if "FAIL" not in line:
                outfile.write(line)
            else:
                pass
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                    #Create BAM index using SAMTOOLS
def bam_index(bam):
    os.system(samtools+" index "+bam)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                    #Create MD5 Sums
def md5calc(input):
    print("Generating MD5 Sums for the results.\nIt is recommended that you copy this file and keep it somewhere safe for any future reference.")
    os.system(md5sum+" "+input+" > "+input+"_MD5_SAVE-ME.txt")
    files.append(input+"_MD5_SAVE-ME.txt")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                    # BAM to EFF
def vars_ann_bam(bam):
    if os.path.isfile(bam):
        print("File seems valid.\nProceeding")
        print("Generating Quality report...")
        os.system(fastqc+" -t "+str(threads)+" -f bam_mapped -extract "+bam)
        if os.path.isfile(bam + ".bai"):
            print("\nLooks like we are all set for the file '{}'\nLet's get to work!\nSummoning GATK HaplotypeCaller!\nRelax! Have a mug of coffee.\nThis might take some time.".format(bam))
            os.system(gatk + " HaplotypeCaller --native-pair-hmm-threads "+str(threads) +" -R " + genome_assembly + " -I " + bam + " -D " + dbsnp + " -L " + bed_file + " -O " + bam[:-4] + ".vcf")
            vcf = bam[:-4] + ".vcf"
            md5calc(vcf)
            if os.path.isfile(vcf):
                print("GATK HaplotypeCaller has generated the VCF file:\n{}\nSummoning GATK SelectVariants".format(vcf))
                os.system(gatk + " SelectVariants -R " + genome_assembly + " -V " + vcf + " --select-type SNP -O " + vcf[:-4] + "_SNVs.vcf")
                os.system(gatk + " SelectVariants -R " + genome_assembly + " -V " + vcf + " --select-type INDEL -O " + vcf[:-4] + "_INDELs.vcf")
                print("\n\nVariants split.\nLet's filter them based on the provided hard filters.\nSummoning GATK VariantFiltration!")
                os.system(gatk + " VariantFiltration -R " + genome_assembly + " -V " + vcf[:-4] + "_SNVs.vcf " + filter_snp + " -O " + vcf[:-4] + "_SNVs_filtered.vcf")
                os.system(gatk + " VariantFiltration -R " + genome_assembly + " -V " + vcf[:-4] + "_INDELs.vcf " + filter_indel + " -O " + vcf[:-4] + "_INDELs_filtered.vcf")
                only_pass(vcf[:-4]+"_SNVs_filtered.vcf")
                only_pass(vcf[:-4]+"_INDELs_filtered.vcf")
                os.system(gatk+" MergeVcfs -I "+vcf[:-4]+"_SNVs_filtered_PASS.vcf"+" -I "+vcf[:-4]+"_INDELs_filtered_PASS.vcf"+" -O "+vcf[:-4]+"_ALL_PASS.vcf")
                md5calc(vcf[:-4]+"_ALL_PASS.vcf")
                # Take the variants which passed the filter.
                os.system(effector + vcf[:-4]+"_ALL_PASS.vcf > " + vcf[:-4] + "_ALL_PASS_ANN.vcf")
                md5calc(vcf[:-4] + "_ALL_PASS_ANN.vcf")
            else:
                print("GATK HaplotypeCaller Failed to create a VCF file!\nCheck the GATK Error above and try to fix it.\nSkipping:\n{}\n--------\n.".format(bam))
                pass
        else:
            print("Your BAM file '{}' lacks an index!\nDon't worry.\nI'll handle it.\nSummoning SAMTOOLS!\n".format(bam))
            bam_index(bam)
            print("\nLooks like we are all set for the file '{}'\nLet's get to work!\nSummoning GATK HaplotypeCaller!\nRelax! Have a mug of coffee.\nThis might take some time.".format(
                    bam))
            os.system(
                gatk + " HaplotypeCaller --native-pair-hmm-threads "+str(threads) +" -R " + genome_assembly + " -I " + bam + " -D " + dbsnp + " -L " + bed_file + " -O " + bam[
                                                                                                                               :-4] + ".vcf")
            vcf = bam[:-4] + ".vcf"
            md5calc(vcf)
            if os.path.isfile(vcf):
                print("GATK HaplotypeCaller has generated the VCF file:\n{}\nHold your horses!\nWe need to filter it first.\nBut first, let me separate SNVs and INDELs\nSummoning GATK SelectVariants")
                os.system(gatk + " SelectVariants -R " + genome_assembly + " -V " + vcf + " --select-type SNP -O " + vcf[:-4] + "_SNVs.vcf")
                os.system(gatk + " SelectVariants -R " + genome_assembly + " -V " + vcf + " --select-type INDEL -O " + vcf[:-4] + "_INDELs.vcf")
                print("\n\nVariants split.\nLet's filter them based on the prvided hard filters.\nSummoning GATK VariantFiltration!")
                os.system(gatk + " VariantFiltration -R " + genome_assembly + " -V " + vcf[:-4] + "_SNVs.vcf " + filter_snp + " -O " + vcf[:-4] + "_SNVs_filtered.vcf")
                os.system(gatk + " VariantFiltration -R " + genome_assembly + " -V " + vcf[:-4] + "_INDELs.vcf " + filter_indel + " -O " + vcf[:-4] + "_INDELs_filtered.vcf")
                only_pass(vcf[:-4]+"_SNVs_filtered.vcf")
                only_pass(vcf[:-4]+"_INDELs_filtered.vcf")
                os.system(gatk+" MergeVcfs -I "+vcf[:-4]+"_SNVs_filtered_PASS.vcf"+" -I "+vcf[:-4]+"_INDELs_filtered_PASS.vcf"+" -O "+vcf[:-4]+"_ALL_PASS.vcf")
                md5calc(vcf[:-4]+"_ALL_PASS.vcf")
                os.system(effector + vcf[:-4]+"_ALL_PASS.vcf > " + vcf[:-4] + "_ALL_PASS_ANN.vcf")
                md5calc(vcf[:-4] + "_ALL_PASS_ANN.vcf")
            else:
                print("GATK HaplotypeCaller Failed to create a VCF file!\nCheck the GATK Error above and try to fix it.\nSkipping:\n{}\n--------\n".format(bam))
                pass
    else:
        print("\nInvalid file '{}' given.\nSkipping...\nDon't try to fool me lad!\n".format(bam))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
                    #FASTQ to EFF
def vars_ann_fq(fastq):
    if fastq.endswith('.fastq'):
        if "_R2" in fastq:
            print("Generating Quality report...")
            os.system(fastqc+" -t "+str(threads)+" -f fastq -extract "+fastq)            
            print("It's the second file for the paired read.\nIt'll be run alongwith its R1 counterpart.")
            pass
        elif "_R1" in fastq:
            print("Generating Quality report...")
            os.system(fastqc+" -t "+str(threads)+" -f fastq -extract "+fastq)            
            if os.path.isfile(fastq[:-7]+"2.fastq"):
                os.system(bwa+" mem -t "+str(threads)+" -R "+readgroup+" "+genome_assembly+" "+fastq+" "+fastq[:-7]+"2.fastq > "+fastq[:-8]+"align.sam")
                os.system(picard+" SortSam SORT_ORDER=coordinate I="+fastq[:-8]+"align.sam O="+fastq[:-8]+"align.bam")
                print("Deleting SAM file....")
                os.system("rm "+fastq[:-8]+"align.sam")
                md5calc(fastq[:-8]+"align.bam")
                vars_ann_bam(fastq[:-8]+"align.bam")
            else:
                print("The corresponding file: {}, not found.\nAre you sure it's there?\nCheck.\nx-x-x-x Skipping this file for now x-x-x-x".format(fastq[:-7]+"2.fastq"))
                pass
        else:
            print("Generating Quality report...")
            os.system(fastqc+" -t "+str(threads)+" -f fastq -extract "+fastq)            
            os.system(bwa+" mem -t "+str(threads)+" -R "+readgroup+" "+genome_assembly+" "+fastq+" > "+fastq[:-6]+"align.sam")
            os.system(picard+" SortSam SORT_ORDER=coordinate I="+fastq[:-6]+"align.sam O="+fastq[:-6]+"align.bam")
            print("Deleting SAM file....")
            os.system("rm "+fastq[:-6]+"align.sam")            
            md5calc(fastq[:-6]+"align.bam")
            vars_ann_bam(fastq[:-6]+"align.bam")
    elif fastq.endswith('.fq'):
        if "_R2" in fastq:
            print("Generating Quality report...")
            os.system(fastqc+" -t "+str(threads)+" -f fastq -extract "+fastq)
            print("It's the second file for the paired read.\nIt'll be run alongwith its R1 counterpart.")
            pass
        elif "_R1" in fastq:
            print("Generating Quality report...")
            os.system(fastqc+" -t "+str(threads)+" -f fastq -extract "+fastq)            
            if os.path.isfile(fastq[:-4] + "2.fastq"):
                os.system(bwa+" mem -t "+str(threads)+" -R "+readgroup+" "+genome_assembly+" "+fastq+" "+fastq[:-4]+"2.fastq > "+fastq[:-5]+"_align.sam")
                os.system(picard+" SortSam SORT_ORDER=coordinate I="+fastq[:-5]+"_align.sam O="+fastq[:-5]+"_align.bam")
                print("Deleting SAM file....")
                os.system("rm "+fastq[:-5]+"align.sam")                 
                md5calc(fastq[:-5]+"_align.bam")
                vars_ann_bam(fastq[:-5]+"_align.bam")
            else:
                print("The corresponding file: {}, not found.\nAre you sure it's there?\nCheck.\nx-x-x-x Skipping this file for now x-x-x-x".format(fastq[:-4]+"2.fastq"))
                pass
        else:
            print("Generating Quality report...")
            os.system(fastqc+" -t "+str(threads)+" -f fastq -extract "+fastq)            
            os.system(bwa+" mem -t "+str(threads)+" -R "+readgroup+" "+genome_assembly+" "+fastq+" > "+fastq[:-3]+"_align.sam")
            os.system(picard+" SortSam SORT_ORDER=coordinate I="+fastq[:-3]+"_align.sam O="+fastq[:-3]+"_align.bam")
            print("Deleting SAM file....")
            os.system("rm "+fastq[:-3]+"align.sam")            
            md5calc(fastq[:-3]+"_align.bam")
            vars_ann_bam(fastq[:-3]+"_align.bam")
    else:
        print("The file {}, doesn't end in FASTQ/FQ.\nIf it is a FASTQ file, kindly check it's extension.\nSkipping the file now.".format(fastq))
        pass
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#==============================================Multiplexing=====================================================================================================================================#

total_files = ((len(sys.argv)))
if total_files-1 == 0:
    print("You've given me 0 files!\nWhat do you want me to do with these many files ?\nExiting...")
    exit()
else:
    print("You've given me {} file(s) to process.\nProceeding...".format(total_files-1))
count = 1
#----------Don't ask. Just for fun.------------#
refreshment = "a mug of coffee"
if total_files >0 and total_files <=2:
    refreshment = "a mug of coffee"
elif total_files >=3 and total_files <=5:
    refreshment = "a walk, and a coffee,"
elif total_files >=6 and total_files <=10:
    refreshment = "a quick nap"
elif total_files >10:
    refreshment= "a day break"
print("Have {} while I do the work for you.".format(refreshment))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
file=sys.argv[1]
while count < total_files:
    print("================================\n{0}. {1}\n================================".format(count,sys.argv[count]))
    if (sys.argv[count]).endswith('.bam'):
        print("It's BAM",(sys.argv[count]))
        vars_ann_bam(sys.argv[count])
        count += 1
    elif (sys.argv[count]).endswith('.fastq') or (sys.argv[count]).endswith('.fq'):
        print("It's FASTQ",(sys.argv[count]))
        vars_ann_fq(sys.argv[count])
        count+=1
    elif (sys.argv[count]).endswith('.fastq.gz') or (sys.argv[count]).endswith('.fq.gz'):
        file = (sys.argv[count])[:-3]
        print("Compressed FASTQ file {} provided.\nI will decompress it in place and DELETE the compressed one to save disk space.\nProceeding...".format(sys.argv[count]))
        if "_R1.fastq" not in sys.argv[count] and "_R1.fq" not in sys.argv[count] and "_R2" not in sys.argv[count]:
            os.system(gunzip+" -v "+sys.argv[count])
        elif "_R1.fq" in sys.argv[count]:
            os.system(gunzip + " -v " + sys.argv[count])
            os.system((gunzip + " -v " + sys.argv[count])[:-7]+"2.fq.gz")
        elif "_R1.fastq" in sys.argv[count]:
            os.system(gunzip + " -v " + sys.argv[count])
            os.system((gunzip + " -v " + sys.argv[count])[:-10]+"2.fastq.gz")
        elif "_R2" in sys.argv[count]:
            print("R2 file for R1 given.\nI must have already processed that file if you gave me R1.\nSkipping this file..\n")
            pass
        if os.path.isfile(file):
            vars_ann_fq(file)
        else:
            print("Uhoh! Well, this is embarassing!\nThe file decompression failed!\nPlease check the error log above.\nSkipping file.")
            pass
        count+=1
    elif (sys.argv[count]).endswith('.bam.gz'):
        print("Compressed BAM file {} provided.\nI will decompress it in place and DELETE the compressed one to save disk space.\nProceeding...".format(sys.argv[count]))
        os.system(gunzip+" -v "+sys.argv[count])
        if os.path.isfile(file):
            vars_ann_bam(file)
        else:
            print("Uhoh! Well, this is embarassing!\nThe file decompression failed!\nPlease check the error log above.\nSkipping file.")
            pass
        count+=1
    else:
        print("The file isn't a FASTQ or BAM or compressed FASTQ or compressed BAM!.\nWasn't I clear enough?\nSkipping the file {}".format(sys.argv[count]))
        count+=1
        pass
#==============================================Closing=====================================================================================================================================#

print(
    "Done!\n===============================\nThank you for using {0}.\nFor any suggestions or complaints, contact:\nAditya Singh\naditya.onco@gmail.com\n----------\nIf I was useful, plesae cite:\n{1}===============================".format(software,citations))
