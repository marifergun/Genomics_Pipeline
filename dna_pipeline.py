import os
import glob
import pandas as pd
from datetime import datetime
import json
import shutil


#specify the paths first

picard_path = "/media/data/AMBRY/picard.jar"
gatk_path = "/media/data/AMBRY/GenomeAnalysisTK.jar"
ref_dir = "/media/data/ref_genome_indexes/"
varscan_path = "/media/data/AMBRY/VarScan.v2.3.9.jar"
dbsnp = "/media/data/ref_genome_indexes/hg19_bundle/dbsnp_138.hg19.vcf"
cosmic = "/media/data/ref_genome_indexes/hg19_bundle/Mills_and_1000G_gold_standard.indels.hg19.vcf"
wd = ""



def system_command_send(command, from_function, th):
    
    logs = {'function': "", 'command': "", 'start_time': "", 'end_time': "", 'threads': "", 'success': 0}
    logs["function"] = from_function
    logs["command"] = command
    logs["start_time"] = str(datetime.now())
    logs["threads"] = th

    try:
        os.system(command)
        logs["end_time"] = str(datetime.now())
        logs["success"] = 1
        write_logs(logs)

    except:
        logs["end_time"] = str(datetime.now())
        logs["success"] = 0
        write_logs(logs)
        return from_function + " give error with this command -> " + command


def write_logs(log):
    import json
    with open('log_file.txt', 'a') as file:
        file.write(json.dumps(log))
        file.write(",")


class BamPipeline(object):
    # workingdirectory,  map_type, sample_type, library_matching_id, thrds
    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds):
        self.working_directory = working_directory
        wd = working_directory
        log_s = {'function': "", 'command': "", 'start_time': "", 'end_time': "", 'threads': "", 'success': 0}
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = thrds
        write_logs(log_s)
        self.bundle_dir = ref_dir + "hg19_bundle"

    def get_fastq(self):
        os.chdir(self.working_directory)
        all_fastq_files = glob.glob("*fastq.gz")

        split_names_v = [os.path.splitext(os.path.splitext(i)[0])[0] for i in all_fastq_files]

        return split_names_v

    def get_info(self, fastq_list):
        sample_ID, germline_dna, index_seq, lanes, pairs_r, n_of_seq = (set() for i in range(6))
        if self.sample_type == "Tumor":
            for i in fastq_list:
                sample_ID.add(i.split("_")[0])
                index_seq.add(i.split("_")[1])
                lanes.add(i.split("_")[2])
                pairs_r.add(i.split("_")[3])
                n_of_seq.add(i.split("_")[4])

            list_with_info = {"Sample_ID": list(sample_ID), "Index": list(index_seq), "Lanes": list(lanes),
                              "Pairs": list(pairs_r), "Number_of_seq": list(n_of_seq)}
            return list_with_info
        elif self.sample_type == "Germline":

            for i in fastq_list:
                sample_ID.add(i.split("_")[0])
                germline_dna.add(i.split("_")[1])
                index_seq.add(i.split("_")[2])
                lanes.add(i.split("_")[3])
                pairs_r.add(i.split("_")[4])
                n_of_seq.add(i.split("_")[5])

            list_with_info = {"Sample_ID": list(sample_ID), "Germline": list(germline_dna), "Index": list(index_seq),
                              "Lanes": list(lanes), "Pairs": list(pairs_r), "Number_of_seq": list(n_of_seq)}
            return list_with_info
        else:
            print("raise error and ask again for a valid sample type")

    def mapping(self, fastq_list, info_dict):
        import re
        import gzip

        RG_SM = info_dict["Sample_ID"][0]
        RG_PL = "Illumina"
        RG_LB = self.library_matching_id
        first_fastq_file_dir = self.working_directory + "/" + fastq_list[0] + ".fastq.gz"
        with gzip.open(first_fastq_file_dir) as f:
            first_line = f.readline()

        flowcell_info = str(first_line).split(":")[2]

        for i in info_dict["Lanes"]:
            for k in info_dict["Number_of_seq"]:
                r1 = re.compile(".*" + i + "_R1_" + k)
                read1 = [s + ".fastq.gz" for s in fastq_list if r1.match(s)]

                r2 = re.compile(".*" + i + "_R2_" + k)
                read2 = [s + ".fastq.gz" for s in fastq_list if r2.match(s)]

                RG_ID = flowcell_info + "." + i[-1]
                RG_PU = flowcell_info + "." + info_dict["Index"][0] + "." + i[-1]
                map_bam = ""
                gene_origin = self.map_type + "_" + info_dict["Sample_ID"][0] + "_" + info_dict["Index"][
                    0] + "_" + i + "_" + k + ".bam"

                if self.map_type == "Bwa":
                    add_read_group = ' -R "@RG\\tID:' + RG_ID + '\\tSM:' + RG_SM + '\\tLB:' + RG_LB + '\\tPL:' + \
                                     RG_PL + '\\tPU:' + RG_PU + '" '

                    map_bam = "bwa mem -t " + self.threads + " " + add_read_group + ref_dir + \
                              "Bwa/ucsc.hg19.fasta " + read1[0] + " " + read2[0] + \
                              " | samtools view -@" + self.threads + " -bS - > " + gene_origin
                elif self.map_type == "Bowtie":

                    add_read_group = " --rg-id " + RG_ID + " --rg SM:" + RG_SM + " --rg LB:" + RG_LB + " --rg PL:" + \
                                     RG_PL + " --rg PU:" + RG_PU

                    map_bam = "bowtie2 -p" + self.threads + add_read_group + " -x " + ref_dir + \
                              "Bowtie/hg_19_bowtie2 -1 " + read1[0] + " -2 " + read2[0] + \
                              " | samtools view -@" + self.threads + " -bS - > " + gene_origin
                else:
                    return "Please specify the map type Bwa/Bowtie "

                system_command_send(map_bam, "mapping", self.threads)

                self.convert_sort(gene_origin)

    def convert_sort(self, sort_gene_origin):

        convert_sort = "samtools view -@" + self.threads + " -bS " + sort_gene_origin + " | samtools sort -@" + \
                       self.threads + " -o SortedBAM_" + sort_gene_origin
        system_command_send(convert_sort, "mapping_function;convert_sort_command", self.threads)

    def merge_bams(self, info_dict):

        all_bam_files = glob.glob("SortedBAM*")
        print(all_bam_files)
        # sorted_bam_files = [os.path.splitext(os.path.splitext(i)[0])[0] for i in all_bam_files]
        inputs_list = ""
        for i in all_bam_files:
            inputs_list = inputs_list + "I=" + i + " "

        ouput_name = self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM.bam"

        merge_command = "java -XX:ParallelGCThreads=" + self.threads + \
                        " -jar " + picard_path + " MergeSamFiles " + inputs_list + \
                        " O=" + ouput_name + " USE_THREADING=true"

        system_command_send(merge_command, "merge_bams", self.threads)

    def mark_duplicate(self):

        merged_bam = glob.glob("*_MergedBAM.bam")
        mark_prefix_removed = self.map_type + "_mdup_removed_"

        picardcommand = "java -XX:ParallelGCThreads=" + self.threads + \
                        " -jar " + picard_path + " MarkDuplicates I=" + merged_bam[0] + \
                        " O=" + mark_prefix_removed + "_" + merged_bam[0] + \
                        " M=marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true"
        system_command_send(picardcommand, "mark_duplicate", self.threads)

    def gatk_realign_target_creator(self):

        bamstr = self.map_type + "_mdup_removed*.bam"
        print(bamstr)
        lastbam = glob.glob(bamstr)
        bcal = "java -jar " + gatk_path + " -T RealignerTargetCreator -nt " + \
               self.threads + " -R " + self.bundle_dir + "/ucsc.hg19.fasta -known " + \
               self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf -I " + lastbam[0] + \
               " -o realign_target.intervals"
        system_command_send(bcal, "GATK_RealignTargetCreator", self.threads)

    def gatk_indel_realigner(self):

        bamstr = self.map_type + "_mdup_removed*.bam"
        lastbam = glob.glob(bamstr)

        realigned_last_bam = "IndelRealigned_" + lastbam[0]
        bcal = "java -jar " + gatk_path + " -T IndelRealigner -R " + self.bundle_dir + "/ucsc.hg19.fasta -known " + \
               self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + \
               " -targetIntervals realign_target.intervals --noOriginalAlignmentTags -I " + lastbam[0] + " -o " + \
               realigned_last_bam

        system_command_send(bcal, "GATK_IndelRealigner", self.threads)

    def gatk_base_recalibrator(self):

        bamstr = "IndelRealigned_*.bam"
        lastbam = glob.glob(bamstr)
        basequalityscore = str(lastbam[0]).split(".")[0] + "_bqsr.grp"
        nct = " -nct " + str(self.threads)
        bcal = "java -jar " + gatk_path + nct + " -T BaseRecalibrator -R " + self.bundle_dir + "/ucsc.hg19.fasta -I " +\
               lastbam[0] + " -knownSites " + self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + \
               " -o " + basequalityscore
        system_command_send(bcal, "GATK_BaseRecalibrator", self.threads)

    def gatk_print_reads(self):

        bamstr = "IndelRealigned_*.bam"
        lastbam = glob.glob(bamstr)
        bqsr = glob.glob("*.grp")[0]
        nct = " -nct " + str(self.threads)
        aftercalibratorBam = "Completeted_BaseCalibrator_" + lastbam[0]
        bcal = "java -jar " + gatk_path + nct + " -T PrintReads -R " + self.bundle_dir + "/ucsc.hg19.fasta -I " + \
               lastbam[0] + " --BQSR " + bqsr + " -o " + aftercalibratorBam

        system_command_send(bcal, "GATK_PrintReads", self.threads)

    def run_gatks(self):
        self.gatk_realign_target_creator()
        self.gatk_indel_realigner()
        self.gatk_base_recalibrator()
        self.gatk_print_reads()

    def run_pipeline(self):

        fastqs = self.get_fastq()
        print(fastqs)
        info = self.get_info(fastqs)
        print(info)
        self.mapping(fastqs, info)
        self.merge_bams(info)
        self.mark_duplicate()
        self.run_gatks()
        self.create_folder()

        return True

    def create_folder(self):

        mk_dir = self.working_directory + "/" + self.map_type
        os.mkdir(mk_dir)
        all_files = glob.glob("*.*")
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)


class VariantCall(object):

    def __init__(self, variant_caller, thrds, map_type, germline_bam, germline_realign):

        self.working_directory = wd + "/" + map_type
        os.chdir(self.working_directory)
        self.v_caller = variant_caller
        self.threads = thrds
        self.map_type = map_type
        self.ref_dir = ref_dir + "/hg19_bundle/ucsc.hg19.fasta"
        normal_bam = glob.glob("Completeted_BaseCalibrator_*.bam")
        normal_realign = glob.glob("realign_target.intervals")
        self.normal_bam = normal_bam[0]
        self.germline_bam = germline_bam
        self.normal_realign = normal_realign[0]
        self.germline_realign = germline_realign
        self.realign_target = self.normal_realign + " " + self.germline_realign
        self.output_vcf = variant_caller + "_ouput.vcf"

    def run_pipeline(self):
        if self.v_caller == "Mutect2":
            self.mutect_caller()
        elif self.v_caller == "Varscan":
            self.varscan_caller()
        return True

    def mutect_caller(self):
        nct = " -nct " + self.threads
        command = "java -jar " + gatk_path + " -T MuTect2 " + nct + " -R " + self.ref_dir + " -I:tumor " + \
                  self.germline_bam + " -I:normal " + self.normal_bam + " --dbsnp " + dbsnp + " --cosmic " + \
                  cosmic + " -L " + self.normal_realign + " -L " + self.germline_realign + " -o " + self.output_vcf
        system_command_send(command, "mutect_caller", self.threads)

    def varscan_caller_step1(self):
        command = "samtools mpileup -f " + self.ref_dir + " -q 1 -B " + self.normal_bam + " " + \
                  self.germline_bam + " > intermediate_mpileup.pileup"

        system_command_send(command, "varscan_caller_step1", self.threads)
        intermediate_mpileup = glob.glob("intermediate_mpileup.pileup")
        print(command)
        return intermediate_mpileup

    def varscan_caller_step2(self, intermediate_mpileup):
        cwd = os.getcwd()
        cwd += "/output.basename"
        command = "java -jar " + varscan_path + " somatic " + intermediate_mpileup + " " + cwd + \
                  " --mpileup 1 --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6 --min-var-freq 0.10 " +\
        "--min-freq-for-hom 0.75 --normal-purity 1.0 --tumor-purity 1.00 --p-value 0.99 --somatic-p-value 0.05 " + \
                  "--strand-filter 0 --output-vcf"

        system_command_send(command, "varscan_caller_step2", self.threads)
        intermediate_varscan_somatic = glob.glob("output.basename*")
        print(command)
        return intermediate_varscan_somatic

    def varscan_caller_step3(self, intermediate_varscan_somatic):
        print(intermediate_varscan_somatic)
        for somatic in intermediate_varscan_somatic:
            command = "java -jar " + varscan_path + " processSomatic " + somatic + " --min-tumor-freq 0.10 " + \
                      "--max-normal-freq 0.05 --p-value 0.07"
            system_command_send(command, "varscan_caller_step3", self.threads)
        return glob.glob("output.basename*")

    def varscan_caller(self):
        step1 = self.varscan_caller_step1()
        step2 = self.varscan_caller_step2(step1)
        step3 = self.varscan_caller_step3(step2)
        print(step3)


class ReadLogs(object):
    def __init__(self, log_file):
        self.log_file = log_file

    def convert_date(self, dt):
        dtt = dt.split(" ")
        dtt_d = dtt[0].split("-")
        dtt_h = dtt[1].split(":")
        da = datetime(year=int(dtt_d[0]), month=int(dtt_d[1]), day=int(dtt_d[2]), hour=int(dtt_h[0]), minute=int(dtt_h[1]),
                      second=int(float(dtt_h[2])))
        return da

    def read_log_file(self):
        open_log = open(self.log_file, "r")
        logs = open_log.read()
        logs = logs[1:-2]
        aa = logs.split("},{")
        bb = [json.loads("{" + a + "}") for a in aa]
        df = pd.DataFrame.from_dict(bb)
        a = self.convert_date(df.start_time.values.min())
        b = self.convert_date(df.end_time.values.max())
        c = b - a
        d = c.total_seconds()
        k = divmod(d, 3600)[0]
        return df, k