# inputs:
# working_dir - the path for the files
# final_flag - a boolean variable, indicates whether we can read the final association file
# assoc_flag - same as the previous but for the association file before filtering for minimal number of associations
# promiscuity_flag - same as before but before filtering for barcode promiscuity
# final_file_name - optional, can be set to be different from the default
# all files must include two columns - bc_seq (barcode sequence) and oligo (the tested CRE)
# output: saves the selected files in the working directory, these files are the input for the QC pipeline


import polars as pl
import sys
print(sys.argv)
working_dir=sys.argv[1]
final_flag = int(sys.argv[2])
assoc_flag = int(sys.argv[3])
promiscuity_flag = int(sys.argv[4])

if len(sys.argv)==6:
    final_file_name = sys.argv[5]

else:
    final_file_name = 'after_associations_filter.csv'

    
if final_flag:
    print(f"{working_dir}/{final_file_name}")
    if final_file_name.lower().endswith(".txt"):
        filtered_df=pl.read_csv(f"{working_dir}/{final_file_name}",separator="\t")
    else:
        filtered_df=pl.read_csv(f"{working_dir}/{final_file_name}")
    filtered_df=filtered_df[['bc_seq','oligo']]
    filtered_df=filtered_df.group_by(['bc_seq','oligo']).len()
    filtered_df=filtered_df.rename({'len':'match_count','bc_seq':'barcode','oligo':'cCRE'})
    filtered_df=filtered_df.to_pandas()
    filtered_df.to_csv(f"/{working_dir}/final_df.csv")
    
if assoc_flag:
    unfiltered_assoc=pl.read_csv(f"{working_dir}/after_promiscuity_filter.csv")
    unfiltered_assoc=unfiltered_assoc[['bc_seq','oligo']]
    unfiltered_assoc=unfiltered_assoc.group_by(['bc_seq','oligo']).len()
    unfiltered_assoc=unfiltered_assoc.rename({'len':'match_count','bc_seq':'barcode','oligo':'cCRE'})
    unfiltered_assoc=unfiltered_assoc.to_pandas()
    unfiltered_assoc.to_csv(f"/{working_dir}/before_min_assoc_df.csv")
    
if promiscuity_flag:
    unfiltered_prom=pl.read_csv(f"{working_dir}/after_barcode_filter.csv")
    unfiltered_prom=unfiltered_prom[['bc_seq','oligo']]
    unfiltered_prom=unfiltered_prom.group_by(['bc_seq','oligo']).len()
    unfiltered_prom=unfiltered_prom.rename({'len':'match_count','bc_seq':'barcode','oligo':'cCRE'})
    unfiltered_prom=unfiltered_prom.to_pandas()
    unfiltered_prom.to_csv(f"/{working_dir}/before_prom_df.csv")