import pandas as pd
from datetime import date
today = date.today()
# import pprint
# import os
# try:
#   path = "G:/NGS_Service/SEQ RUNS/2024-06-24-NextSeq"
# except: 
#   path = input("Please provide the path of the source directory: ")
# os.chdir(path)

# Create dictionaries for the data
ds = []
ss1 = [] # ss = [], not needed, only kept ss1 name for legacy reasons
mq = []
mq_header = ["SampleID","Filename","File_type","Encoding","Total_Sequences","Total_Bases_Mbp","Sequences_flagged_as_poor_quality","Sequence_length","GC","total_deduplicated_percentage","avg_sequence_length","median_sequence_length","basic_statistics","per_base_sequence_quality"," per_tile_sequence_quality","per_sequence_quality_scores","per_base_sequence_content","per_sequence_gc_content","per_base_n_content","sequence_length_distribution","sequence_duplication_levels","overrepresented_sequences","adapter_content"]
# mq1 = []

# Append all columns, prepare ds
with open("Demultiplex_Stats.csv") as ds_source:
    header = ds_source.readline().strip().split(",")
    for row in ds_source:
        data = row.strip().split(",")
        for i in range(len(data)):
            sample = {header[0] : data[0], header[1] : data[1], header[2] : data[2], header[3] : data[3], header[4] : data[4], header[5] : data[5], header[6] : data[6], header[7] : data[7], header[8] : data[8], header[9] : data[9], header[10] : data[10]}
        ds.append(sample)
# return and delete value of undetermined reads
undetermined = ds.pop(-1)

# Append PASS/FAILED to a yield_seq list
threshold = 75000
for i in range(len(ds)):
    if int(ds[i]["# Reads"]) > threshold:
        ds[i]["Yield"] = "PASS"
    elif int(ds[i]["# Reads"]) < threshold:
        ds[i]["Yield"] = ("FAILED")
    elif ds[i]["# Reads"] == "NA":
        ds[i]["Yield"] = ("NA")

# Append only selected columns, prepare ss1
with open("QC_SampleSheet.csv") as ss_source:
    # remove first row
    ss_source.readline()
    # read header
    header = ss_source.readline().strip().split(",")
    # read data
    for row in ss_source:
        data = row.strip().split(",")
        for i in range(len(data)):
            sample = {"SampleID" : data[0], header[1] : data[1], header[2] : data[2], header[3] : data[3], header[4] : data[4]}
        ss1.append(sample)

# Function to rename samples
def rename_sample(data_col_num):
        data_x = data[data_col_num].split("_")
        name_to_keep = len(data_x)-4
        for i in range(len(data_x)):
            while len(data_x) > name_to_keep:
                data_x.pop(-1)
        data[data_col_num] = "_".join(data_x)

# Append only selected columns, strip data to format, prepare mq
with open("QC_multiqc_fastqc.txt") as mq_source:
    for row in mq_source:
        data = row.strip().split("\t")
        rename_sample(0)
        rename_sample(1)
        data[1] = data[1]+".fastq.gz"
        data[5] = data[5].strip(" Mbp")
        if "kbp" in data[5]:
            data[5] = data[5].strip(" kbp")
            data[5] = float(data[5])/1000
        sample = {mq_header[0] : data[0], mq_header[1] : data[1], mq_header[2] : data[2], mq_header[3] : data[3], mq_header[4] : data[4], mq_header[5] : data[5], mq_header[6] : data[6], mq_header[7] : data[7], mq_header[8] : data[8], mq_header[9] : data[9], mq_header[10] : data[10], mq_header[11] : data[11], mq_header[12] : data[12], mq_header[13] : data[13], mq_header[14] : data[14], mq_header[15] : data[15], mq_header[16] : data[16], mq_header[17] : data[17], mq_header[18] : data[18], mq_header[19] : data[19], mq_header[20] : data[20], mq_header[21] : data[21], mq_header[22] : data[22]}        
        mq.append(sample)

# Prepare mqb from mq
# rename GC to GC mean
mqb = []
for i in range(len(mq)-1):
    if mq[i]["SampleID"] == mq[i+1]["SampleID"]:
        merged_sample = {mq_header[0] : mq[i+1]["SampleID"], mq_header[4] : float(mq[i+1]["Total_Sequences"])+float(mq[i]["Total_Sequences"]), mq_header[5] : float(mq[i+1]["Total_Bases_Mbp"])+float(mq[i]["Total_Bases_Mbp"]), mq_header[8] : (float(mq[i+1]["GC"])+float(mq[i]["GC"]))/2, mq_header[13] : (mq[i+1]["per_base_sequence_quality"])+", "+(mq[i]["per_base_sequence_quality"]), mq_header[22] : (mq[i+1]["adapter_content"])+", "+(mq[i]["adapter_content"])}
        # print(merged_sample)
        mqb.append(merged_sample)

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

def xls_to_list():
    # optionally, usecols="A:L" from R script
    sets_to_read = [["A" , 2], ["B", 3], ["C", 3], ["D", 3]]
    for i in range(len(sets_to_read)):
        filename = "Set"+(sets_to_read[i][0])+"_16s_Targeted_Seq_Workflow.xlsx"
        try:
            setx = (pd.read_excel(filename, sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
        except:
            error = filename+" not found. Please provide the actual name (eg. Set"+(sets_to_read[i][0])+".xlsx): "
            file = input(error)
            setx = (pd.read_excel(file, sheet_name = "Bioinfo", usecols="A:W")).values.tolist()
        # print(setx[:5])
        # print(type(setx))
        del setx[:sets_to_read[i][1]]
        set_to_list(setx, sets_to_read[i][0])
xls_to_list()

# Prepare dictionary from sets list
sets_header = sets[0]
# Rename original "Sample" to "SampleID"
sets_header[0] = "SampleID" 
sets_list = []
for i in range(len(sets)-1):
    for j in range(len(sets_header)):
        sets_dict = dict(zip(sets_header, sets[i+1]))    
    sets_list.append(sets_dict)

# Creating dataframes
ss1_df = pd.DataFrame(ss1)
mqb_df = pd.DataFrame(mqb)
mq_df = pd.DataFrame(mq)
sets_list_df = pd.DataFrame(sets_list)
ds_df = pd.DataFrame(ds)

# Merging dataframes
ds1_df = pd.merge(ss1_df, ds_df, on="SampleID")
qc_df = pd.merge(sets_list_df, ds1_df, on="SampleID") # sets + ds1
qc2_df = (pd.merge(qc_df, mqb_df, on="SampleID")).sort_values("SampleID") # qc + mqb
qc2_df.rename(columns={"GC" : "GC mean"}, inplace=True)
qc3_df = qc2_df.drop(["Total_Sequences", "Total_Bases_Mbp", "GC mean", "per_base_sequence_quality", "adapter_content"], axis=1) 


# Exporting dataframes to csv
multiQC_name = "QC-"+str(today)+"-NextSeq_multiQC.csv"
with open(multiQC_name, "w") as qc_output:
    qc2_df.to_csv(qc_output, index=False, lineterminator="\n")

QC_name = "QC-"+str(today)+"-NextSeq.csv"
with open(QC_name, "w") as qc_output:
    qc_df.to_csv(qc_output, index=False, lineterminator="\n")
    # qc3_df.to_csv(qc_output, index=False, lineterminator="\n")

X_Reads = []
for i in range(len(ds)-1):
    reads_int = int(ds[i+1]["# Reads"])
    X_Reads.append(reads_int)

# Create summary statistics
sum_col = ["min", "max", "passed", "failed", "pass %", "fail %"]
min_reads = min(X_Reads)
max_reads = max(X_Reads)
passed = 0
failed = 0
for i in range(len(ds)):
    if ds[i]["Yield"] == "PASS":
        passed += 1
    elif ds[i]["Yield"] == "FAILED":
        failed += 1

pass_perc = round(passed / (passed + failed) * 100, 3)
fail_perc = round(failed / (passed + failed) * 100, 3)
sum = {sum_col[0]: min_reads, sum_col[1]: max_reads, sum_col[2]: passed, sum_col[3]: failed, sum_col[4]: pass_perc, sum_col[5]: fail_perc}
sum_df = pd.DataFrame([sum])
print(sum)

# Exporting summary statistics to csv
QC_stat_name = "QC-"+str(today)+"-NextSeq_stat.csv"
with open(QC_stat_name, "w") as qc_output:
    sum_df.to_csv(qc_output, index=False, lineterminator="\n")
