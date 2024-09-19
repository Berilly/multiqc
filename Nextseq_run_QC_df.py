import pandas as pd
from datetime import date
today = date.today()
import random
# import os
# try:
#   path = "G:/NGS_Service/SEQ RUNS/"+today+"-NextSeq"
#   os.chdir(path)
# except: 
#   path = input("Please provide the path of the source directory: (Default path = N:/NGS_Service/SEQ RUNS/-today(yyyy-mm-dd)-NextSeq)")
#   os.chdir(path)

# Create dataframes from the data
ds_df = pd.read_csv("Demultiplex_Stats.csv")
undetermined_reads = str(ds_df.loc[384, "# Reads"])
ds_df = ds_df.drop(index=384)

# Append PASS/FAILED to ds_df
threshold = 75000
ds_df["Yield"] = ds_df["# Reads"].apply(lambda x: "PASS" if int(x) > threshold else "FAILED" if int(x) < threshold else "NA")

# Create dataframes from the data
ss_df = pd.read_csv("QC_SampleSheet.csv",  skiprows=1, usecols=["Sample_ID", "ProjectName"])
ss_df.rename(columns={"Sample_ID" : "SampleID"}, inplace=True)

# Create a header for the multiqc data
mq_header = ["SampleID","Filename","File_type","Encoding","Total_Sequences","Total_Bases_Mbp","Sequences_flagged_as_poor_quality","Sequence_length","GC","total_deduplicated_percentage","avg_sequence_length","median_sequence_length","basic_statistics","per_base_sequence_quality"," per_tile_sequence_quality","per_sequence_quality_scores","per_base_sequence_content","per_sequence_gc_content","per_base_n_content","sequence_length_distribution","sequence_duplication_levels","overrepresented_sequences","adapter_content"]
mq_df = pd.read_csv("QC_multiqc_fastqc.txt", sep="\t", names=mq_header)
for i in range(len(mq_df["Total_Bases_Mbp"])):
    mq_df.loc[i, "Total_Bases_Mbp"] = mq_df.loc[i, "Total_Bases_Mbp"].strip(" Mbp")
    if "kbp" in mq_df.loc[i, "Total_Bases_Mbp"]:
        mq_df.loc[i, "Total_Bases_Mbp"] = mq_df.loc[i,"Total_Bases_Mbp"].strip(" kbp")
        mq_df.loc[i, "Total_Bases_Mbp"] = float(mq_df.loc[i, "Total_Bases_Mbp"])/1000
mq_df_renamed = mq_df.replace(" Mbp", "", regex=True)
mq_df_renamed = mq_df_renamed.replace(" kbp", "", regex=True)
mq_df_renamed["Total_Bases_Mbp"] = mq_df_renamed["Total_Bases_Mbp"].astype(float)
mq_df_renamed["Total_Sequences"] = mq_df_renamed["Total_Sequences"].astype(int)
mq_df_renamed["GC"] = mq_df_renamed["GC"].astype(float)

# Function to rename samples
def rename_sample(df, data_col_name):
        for i in range(len(df[data_col_name])):
            data_x = df[data_col_name][i].split("_")
            data_x = data_x[:len(data_x)-4]
            data_x = "_".join(data_x)
            df.replace(df.loc[i,data_col_name], data_x, inplace=True)
renamed = rename_sample(mq_df_renamed, "SampleID")

# Function to merge samples in mq
mqb_df_temp = mq_df_renamed[["SampleID", "Total_Sequences", "Total_Bases_Mbp", "GC", "per_base_sequence_quality","adapter_content"]].copy()
for entry in mqb_df_temp:
     mqb_df = mqb_df_temp.groupby("SampleID").agg({"Total_Sequences" : "sum", "Total_Bases_Mbp" : "sum", "GC" : "mean", "per_base_sequence_quality" : list, "adapter_content" : list})
mqb_df.reindex()
mqb_df.rename(columns={"GC" : "GC mean"}, inplace=True)
mqb_df["per_base_sequence_quality"] = mqb_df["per_base_sequence_quality"].apply(lambda x: str(x).strip("[]").replace("'", ""))
mqb_df["adapter_content"] = mqb_df["adapter_content"].apply(lambda x: str(x).strip("[]").replace("'", ""))
# print(mqb_df.head(5))

# import and merge of excel lists
def sets_import_to():
    sets_df_temp = pd.DataFrame()
    # optionally, usecols="A:L" from R script
    sets_to_read = ["A", "B", "C", "D"]
    for i in range(len(sets_to_read)):
        setname = "Set"+sets_to_read[i]
        filename = "Set"+(sets_to_read[i])+"_16s_Targeted_Seq_Workflow.xlsx"
        try:
            setx_df = (pd.read_excel(filename, sheet_name = "Bioinfo", usecols="A:W", skiprows=3, converters={"Library Plate": str}))
        except:
            error = filename+" not found. Please provide the actual name (eg. Set"+(sets_to_read[i])+"_16s_Targeted_Seq_Workflow.xlsx): "
            file = input(error)
            setx_df = (pd.read_excel(file, sheet_name = "Bioinfo", usecols="A:W", skpirows=3, converters={"Library Plate": str}))
        for j in range(len(setx_df["Library Plate"])):
            setx_df.loc[j, "Library Plate"] = str(setname)
        if len(sets_df_temp) == None:
                sets_df_temp = setx_df
        else:
            sets_df_temp= pd.concat([sets_df_temp, setx_df], axis=0, ignore_index=True)
        sets_df_temp["QuantIT DNA concentration (ng/µl)"] = sets_df_temp["QuantIT DNA concentration (ng/µl)"].apply(lambda x: round(float(x), 4) if x != "QuantIT DNA concentration (ng/µl)" else x)
        sets_df_temp["nM concentration"] = sets_df_temp["nM concentration"].apply(lambda x: round(float(x), 4) if x != "nM concentration" else x)
        # bad programming practice, but atm Conrol and Control info are not used
        sets_df_temp["Control"] = sets_df_temp["Control"].apply(lambda x: "NA" if x != "" else x)
        sets_df_temp["Control info"] = sets_df_temp["Control info"].apply(lambda x: "NA" if x != "" else x)
    return sets_df_temp
sets_df = sets_import_to()
sets_df.rename(columns={"Sample" : "SampleID"}, inplace=True)        
sets_df.rename(columns={"Unnamed: 16" : "...17"}, inplace=True)
sets_df.rename(columns={"Unnamed: 18" : "...19"}, inplace=True)

# Merging dataframes
ds1_df = pd.merge(ss_df, ds_df, on="SampleID")
qc_df = pd.merge(sets_df, ds1_df, on="SampleID") # sets + ds1
qc2_df = (pd.merge(qc_df, mqb_df, on="SampleID")).sort_values("SampleID") # qc + mqb
qc2_df.rename(columns={"# Reads" : "X..Reads"}, inplace=True)
qc2_df.rename(columns={"# Perfect Index Reads" : "X..Perfect.Index.Reads"}, inplace=True)
qc2_df.rename(columns={"# One Mismatch Index Reads" : "X..One.Mismatch.Index.Reads"}, inplace=True)
qc2_df.rename(columns={"% Reads" : "X..Reads.1"}, inplace=True)
qc2_df.rename(columns={"% Perfect Index Reads" : "X..Perfect.Index.Reads.1"}, inplace=True)
qc2_df.rename(columns={"% One Mismatch Index Reads" : "X..One.Mismatch.Index.Reads.1"}, inplace=True)
qc2_df.rename(columns={"% Two Mismatch Index Reads" : "X..Two.Mismatch.Index.Reads.1"}, inplace=True)
qc3_df = qc2_df.drop(["Total_Sequences", "Total_Bases_Mbp", "GC mean", "per_base_sequence_quality", "adapter_content"], axis=1) 

# Exporting dataframes to csv
multiQC_name = "QC-"+str(today)+"-NextSeq_multiQC.csv"
with open(multiQC_name, "w") as qc_output:
    qc2_df.to_csv(qc_output, index=False, lineterminator="\n")

QC_name = "QC-"+str(today)+"-NextSeq.csv"
with open(QC_name, "w") as qc_output:
    qc_df.to_csv(qc_output, index=False, lineterminator="\n")

X_Reads = []
for i in range(len(ds_df)):
    reads_int = int(ds_df.loc[i, "# Reads"])
    X_Reads.append(reads_int)

# Create summary statistics
sum_col = ["min", "max", "passed", "failed", "pass %", "fail %"]
min_reads = min(X_Reads)
max_reads = max(X_Reads)
passed = int(ds_df["Yield"].value_counts()["PASS"])
failed = int(ds_df["Yield"].value_counts()["FAILED"])

pass_perc = round(passed / (passed + failed) * 100, 3)
fail_perc = round(failed / (passed + failed) * 100, 3)
sum = {sum_col[0]: min_reads, sum_col[1]: max_reads, sum_col[2]: passed, sum_col[3]: failed, sum_col[4]: pass_perc, sum_col[5]: fail_perc, "Undetermined Reads" : undetermined_reads}
disclaimer = "Note: Undetermined reads are not included in the min/max summary statistics."
print(sum)
print(disclaimer)
sum_df = pd.DataFrame([sum])
print(sum_df)

# Exporting summary statistics to csv
QC_stat_name = "QC-"+str(today)+"-NextSeq_stat.csv"
with open(QC_stat_name, "w") as qc_output:
    sum_df.to_csv(qc_output, index=False, lineterminator="\n")

# Exporting qc2_df to json
QC_json_name = "QC-"+str(today)+"-NextSeq.json"
with open(QC_json_name, "w") as qc_output:
    qc2_df.to_json(qc_output, orient="records", lines=True)