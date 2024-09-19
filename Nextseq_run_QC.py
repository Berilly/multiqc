import pandas as pd
import numpy as np
from datetime import date
today = date.today()
# import os
# try:
#   path = "G:/NGS_Service/SEQ RUNS/2024-06-24-NextSeq"
# except: 
#   path = input("Please provide the path of the source directory: ")
# os.chdir(path)

# Create lists for the data
ds = []
ss1 = [] # ss = [], not needed, only kept ss1 name for legacy reasons
mq = []
mq_header = ["Sample","Filename","File_type","Encoding","Total_Sequences","Total_Bases_Mbp","Sequences_flagged_as_poor_quality","Sequence_length","GC","total_deduplicated_percentage","avg_sequence_length","median_sequence_length","basic_statistics","per_base_sequence_quality"," per_tile_sequence_quality","per_sequence_quality_scores","per_base_sequence_content","per_sequence_gc_content","per_base_n_content","sequence_length_distribution","sequence_duplication_levels","overrepresented_sequences","adapter_content"]
mq.append(mq_header)
mq1 = []
mq1_header = ["Sample","Total_Sequences","Total_Bases_Mbp","GC","per_base_sequence_quality","adapter_content"]
mq1.append(mq1_header)

# Append all columns, prepare ds
with open("Demultiplex_Stats.csv") as ds_source:
    for sor in ds_source:
        data = sor.strip().split(",")
        ds.append(data)

# Append PASS/FAILED to a yield_seq list
threshold = 75000
yield_seq = ["Yield"] 
for i in range(len(ds)-1):
    if int(ds[i+1][3]) > threshold:
        yield_seq.append("PASS")
    elif int(ds[i+1][3]) < threshold:
        yield_seq.append("FAILED")
    elif ds[[i+1][3]] == "NA":
        yield_seq.append("NA")

# Append only selected columns, prepare ss1
with open("QC_SampleSheet.csv") as ss_source:
    for sor in ss_source:
        data = sor.strip().split(",")
        ss1.append(data[:2])
del ss1[0]

# Append only selected columns, strip data to format, prepare mq
with open("QC_multiqc_fastqc.txt") as mq_source:
    for sor in mq_source:
        mq1_temp = []
        mq_all_temp = []
        data = sor.strip().split("\t")
        l = len(data[0])-1
        number_of_underscores = 4
        for i in range(l):
            while number_of_underscores != 0:
                if data[0][(l-i)] == "_":
                    number_of_underscores -= 1        
                    data[0] = data[0][:(l-i)]                    
                break                     
        data[1] = data[1].strip(".fastq.gz")
        k = len(data[1])-1
        number_of_underscores = 5
        for i in range(k):
            while number_of_underscores != 0:
                if data[1][(k-i)] == "_":
                    number_of_underscores -= 1        
                    data[1] = data[1][:(k-i)]                    
                break
        data[1] = data[1]+".fastq.gz"

        data[5] = data[5].strip(" Mbp")
        if "kbp" in data[5]:
            data[5] = data[5].strip(" kbp")
            data[5] = float(data[5])/1000
        mq1_temp.append(data[0])
        mq1_temp.append(data[4])
        mq1_temp.append(data[5])
        mq1_temp.append(data[8])
        mq1_temp.append(data[13])
        mq1_temp.append(data[22])
        mq1.append(mq1_temp)
        mq.append(data)

# Prepare mqb from mq
mqb = []

for i in range(len(mq1)):
    mqb_temp = []
    if mq1[i-1][1] == mq1[i][1]:
        mqb_temp = [mq1[i-1][0], float(mq1[i-1][1])+float(mq1[i][1]), float(mq1[i-1][2])+float(mq1[i][2]), (float(mq1[i-1][3])+float(mq1[i][3]))/2, (mq1[i-1][4])+", "+(mq1[i][4]), (mq1[i-1][5])+", "+(mq1[i][5])]
        mqb.append(mqb_temp)


# import and merge of excel lists
sets = []

def set_to_list(set_name, set_id):
    for i in range(len(set_name)): # make a function?
        set_temp = []
        for j in range(len(set_name[i])):
            set_temp.append(set_name[i][j])
        set_temp[8] = set_id
        if set_temp[11] != "QuantIT DNA concentration (ng/µl)":
            set_temp[11] = round(float(set_temp[11]), 4)
        else:
            set_temp[11] = set_temp[11]
        if set_temp[14] != "nM concentration":
            set_temp[14] = round(float(set_temp[14]), 4)
        else:
            set_temp[14] = set_temp[14]
        sets.append(set_temp)

try:
    seta = (pd.read_excel("SetA_16s_Targeted_Seq_Workflow.xlsx", sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
except:
    file = input("SetA file name does not match expected. Please provide the actual name (eg. SetA.xlsx): ")
    seta = (pd.read_excel(file, sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
# delete only two to keep the header
del seta[:2]
set_to_list(seta,"Set A")

try:
    setb = (pd.read_excel("SetB_16s_Targeted_Seq_Workflow.xlsx", sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
except:
    file = input("SetB file name does not match expected. Please provide the actual name (eg. SetB.xlsx): ")
    setb = (pd.read_excel(file, sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
del setb[:3]
set_to_list(setb, "Set B")

try:
    setc = (pd.read_excel("SetC_16s_Targeted_Seq_Workflow.xlsx", sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
except:
    file = input("SetC file name does not match expected. Please provide the actual name (eg. SetC.xlsx): ")
    setc = (pd.read_excel(file, sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
del setc[:3]
set_to_list(setc, "Set C")

try:
    setd = (pd.read_excel("SetD_16s_Targeted_Seq_Workflow.xlsx", sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
except:
    file = input("SetD file name does not match expected. Please provide the actual name (eg. SetD.xlsx): ")
    setc = (pd.read_excel(file, sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
del setd[:3]
set_to_list(setd, "Set D")

# Creating dataframes
ss1_df = pd.DataFrame(ss1, columns=["Sample_ID","ProjectName"])
mqb_df = pd.DataFrame(mqb, columns=["Sample_ID","Total_Sequences","Total_Bases_Mbp","GC mean","per_base_sequence_quality","adapter_content"])
mq_df = pd.DataFrame(mq, columns=["Sample_ID","Filename","File_type","Encoding","Total_Sequences","Total_Bases_Mbp","Sequences_flagged_as_poor_quality","Sequence_length","GC","total_deduplicated_percentage","avg_sequence_length","median_sequence_length","basic_statistics","per_base_sequence_quality"," per_tile_sequence_quality","per_sequence_quality_scores","per_base_sequence_content","per_sequence_gc_content","per_base_n_content","sequence_length_distribution","sequence_duplication_levels","overrepresented_sequences","adapter_content"])
sets_df = pd.DataFrame(sets, columns=["Sample_ID","Comments","Sample type","Isolation kit","Repetition","Plate","Position","isolation conc.","Library Plate","Library position","PCR1. conc.","QuantIT DNA concentration", "Fragment size", "Purity", "nM concentration", "i7 index", "...17", "i5 index", "...19", "project", "Control", "Control info", "Enzyme used"])
ds_df = pd.DataFrame(ds, columns=["Lane","Sample_ID","Index","Reads","Perfect Index Reads","One Mismatch Index Reads","Two Mismatch Index Reads","% Reads","% Perfect Index Reads","% One Mismatch Index Reads","% Two Mismatch Index Reads"
])
ds_df["Yield"] = yield_seq

# Merging dataframes
ds1_df = pd.merge(ss1_df, ds_df, on="Sample_ID")
qc_df = pd.merge(sets_df, ds1_df, on="Sample_ID") # sets + ds1
qc2_df = (pd.merge(qc_df, mqb_df, on="Sample_ID")).sort_values("Sample_ID") # qc + mqb
qc3_df = qc2_df.drop(["Total_Sequences", "Total_Bases_Mbp", "GC mean", "per_base_sequence_quality", "adapter_content"], axis=1) 

# Exporting dataframes to csv
multiQC_name = "QC-"+str(today)+"-NextSeq_multiQC.csv"
with open(multiQC_name, "w") as qc_output:
    qc2_df.to_csv(qc_output, index=False, lineterminator="\n")

QC_name = "QC-"+str(today)+"-NextSeq.csv"
with open(QC_name, "w") as qc_output:
    qc3_df.to_csv(qc_output, index=False, lineterminator="\n")

yield_t = [np.transpose(yield_seq)]
X_Reads = []
for i in range(len(ds)-1):
    ds[i+1][3] = int(ds[i+1][3])
    X_Reads.append(ds[i+1][3])

# Create summary statistics
sum_col = ["min", "max", "passed", "failed", "pass %", "fail %"]
min_reads = min(X_Reads)
max_reads = max(X_Reads)
passed = 0
failed = 0
for i in range(len(yield_seq)):
    if yield_seq[i] == "PASS":
        passed += 1
    elif yield_seq[i] == "FAILED":
        failed += 1

pass_perc = round(passed / (passed + failed) * 100, 3)
fail_perc = round(failed / (passed + failed) * 100, 3)
sum = {sum_col[0]: min_reads, sum_col[1]: max_reads, sum_col[2]: passed, sum_col[3]: failed, sum_col[4]: pass_perc, sum_col[5]: fail_perc}
sum_df = pd.DataFrame([sum])

# Exporting summary statistics to csv
QC_stat_name = "QC-"+str(today)+"-NextSeq_stat.csv"
with open(QC_stat_name, "w") as qc_output:
    sum_df.to_csv(qc_output, index=False, lineterminator="\n")

print(sum)
