import tkinter as tk
import dna_pipeline as dna
from glob import glob
from os import chdir

class PipelineGui:

    def __init__(self):
        self.parent = tk.Tk()
        self.parent.geometry('650x600')
        self.label0 = tk.Label(self.parent, text="Pipeline GUI", fg="black", font=("Times", 40, "bold italic")) \
            .grid(column=2, row=0, sticky="SE", padx=10)

        self.label1 = tk.Label(self.parent, text="1. Pipeline", fg="black", font=("Times", 30, "bold")) \
            .grid(column=1, row=1, sticky="SE", padx=10)

        self.label2 = tk.Label(self.parent, text="Working Directory: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=2, columnspan=2)
        self.working_directory = tk.Entry(self.parent)
        self.working_directory.grid(column=2, row=2, columnspan=4)

        self.label3 = tk.Label(self.parent, text="Mapping Type: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=3, columnspan=2)

        self.var_maptype = tk.StringVar()

        self.map_type1 = tk.Radiobutton(self.parent, text="Bwa", variable=self.var_maptype, value="Bwa") \
            .grid(column=2, row=3)
        self.map_type2 = tk.Radiobutton(self.parent, text="Bowtie2", variable=self.var_maptype, value="Bowtie") \
            .grid(column=3, row=3)

        self.label4 = tk.Label(self.parent, text="Sample Type: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=4, columnspan=2)

        self.var_sampletype = tk.StringVar()

        self.sample_type1 = tk.Radiobutton(self.parent, text="Tumor", variable=self.var_sampletype, value="Tumor") \
            .grid(column=2, row=4)
        self.sample_type2 = tk.Radiobutton(self.parent, text="Germline", variable=self.var_sampletype, value="Germline") \
            .grid(column=3, row=4)

        self.label5 = tk.Label(self.parent, text="Library ID: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=5, columnspan=2)
        self.library = tk.Entry(self.parent)
        self.library.grid(column=2, row=5, columnspan=4)

        self.label6 = tk.Label(self.parent, text="Threads: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=6, columnspan=2)
        self.threads = tk.Entry(self.parent)
        self.threads.grid(column=2, row=6, columnspan=4)
        self.submit_button = tk.Button(self.parent, text="Next =>", command=self.onSample, width=20).grid(row=7, column=2)
        self.parent.mainloop()

    def OnClick(self):

        mt = self.var_maptype.get()
        st = self.var_sampletype.get()
        wd = self.working_directory.get()
        lb = self.library.get()
        th = self.threads.get()
        vc = ""
        gm = ""
        sample_type_check = False

        if st == "Germline":
            sample_type_check = False
            
        else:
            sample_type_check = True
            vc = self.var_variantcaller.get()
            gm = self.germline.get()

        print(sample_type_check)
        if sample_type_check:
            pipeline1 = dna.BamPipeline(working_directory=wd, map_type=mt, sample_type=st, library_matching_id= lb, thrds=th)
            pipeline1_success = pipeline1.run_pipeline()
            chdir(gm + "/" + mt)
            gm_bam = glob("Completeted_BaseCalibrator_*.bam")
            gm_interval = glob("realign_target.intervals")
            chdir(wd+"/"+mt)
            bam = gm + "/" + mt + "/" + gm_bam[0]
            interval = gm + "/" + mt + "/" + gm_interval[0]
            pipeline2 = dna.VariantCall(variant_caller=vc, thrds=th, map_type=mt, germline_bam=bam, germline_realign=interval)
            pipeline2_success = pipeline2.run_pipeline()
            return pipeline2_success
        else:
            print("------------")
            pipeline1 = dna.BamPipeline(working_directory=wd, map_type=mt, sample_type=st, library_matching_id= lb, thrds=th)
            pipeline1_success = pipeline1.run_pipeline()
            return pipeline1_success

    def onSample(self):

        st = self.var_sampletype.get()
        if st == "Tumor":
            for label in self.parent.grid_slaves():
                if int(label.grid_info()["row"]) == 7:
                    label.grid_forget()
            self.label1_2 = tk.Label(self.parent, text="2. Pipeline", fg="black", font=("Times", 30, "bold")) \
                .grid(column=1, row=7, sticky="SE", padx=10)

            self.var_variantcaller = tk.StringVar()
            self.label7 = tk.Label(self.parent, text="Variant Caller: ", fg="black", font=("Times", 15)) \
                .grid(column=0, row=8, columnspan=2)
            self.variantcaller1 = tk.Radiobutton(self.parent, text="Mutect2", variable=self.var_variantcaller, value="Mutect2") \
                .grid(column=2, row=8)
            self.variantcaller2 = tk.Radiobutton(self.parent, text="Varscan", variable=self.var_variantcaller, value="Varscan") \
                .grid(column=3, row=8)

            self.label8 = tk.Label(self.parent, text="Germline Folder: ", fg="black", font=("Times", 15)) \
                .grid(column=0, row=9, columnspan=2)
            self.germline = tk.Entry(self.parent)
            self.germline.grid(column=2, row=9, columnspan=4)

            self.submit_button_final = tk.Button(self.parent, text="Submit", command=self.OnClick, width=20).grid(row=10, column=2)
        else:
            self.OnClick()


PipelineGui()
