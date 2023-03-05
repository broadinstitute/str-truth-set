# download logs
import asyncio
import collections
from hailtop.config import get_user_config
from hailtop.batch_client.client import BatchClient
import pandas as pd
import re
import time


billing_project = get_user_config().get('batch', 'billing_project', fallback=None)
bc = BatchClient(billing_project=billing_project)

output_rows = []
loci_with_expansion_hunter_errors = []
number_of_data_groups = 8
current_data_group_number = collections.defaultdict(int)

# get batch
for b in list(
        bc.list_batches(limit=1, last_batch_id=7102273)) + list(
        bc.list_batches(limit=1, last_batch_id=7102275)):  # last_batch_id=6958349, 6954002, 6958255, 6958288

    batch_name = b.attributes["name"]
    if not batch_name.startswith("STR Truth Set:"):
        continue
    batch_name_fields = batch_name.split(": ")
    if len(batch_name_fields) != 4:
        print("WARNING: Unexpected batch name: ", batch_name, "Skipping...")
        continue

    print("="*30)
    tool_name = batch_name_fields[1]
    positive_or_negative_loci = batch_name_fields[2]
    input_bam_name = batch_name_fields[3]

    if input_bam_name == "CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram":
        coverage = "40x genome"
    elif "downsampled_to_" in input_bam_name:
        coverage = input_bam_name.split("downsampled_to_")[-1].replace(".bam", "")
    elif input_bam_name.endswith("CHMI_CHMI3_Nex1.cram"):
        coverage = "exome"
    else:
        print(f"WARNING: Unexpected sample: {input_bam_name}. Skipping...")
        continue

    batch_status = b.status()
    print(f"Batch {b.id}:", batch_name, batch_status["time_created"]) #; continue


    for retry in range(3):
        try:
            for job_info in b.jobs():
                if not job_info["name"].startswith("Run "):
                    print(" "*4, "Job " + str(job_info["job_id"]) + ":", job_info["name"], " ---", job_info["state"], " - skipping...")
                    continue
                if job_info["state"] != "Success":
                    print("WARNING: job state is not Success. Skipping job...")
                    continue

                job_obj = b.get_job(job_info["job_id"])
                job_log = job_obj.log()


                #job_input_log_lines = job_log["input"].split("\n")
                print(" "*4, "Job " + str(job_info["job_id"]) + ":", job_info["name"], " ---", job_info["state"])

                if tool_name == "ExpansionHunter":
                    for match in re.finditer("Error on locus spec ([^:]+):.*", job_log["main"]):
                        print(" "*8, match.group(0))
                        loci_with_expansion_hunter_errors.append({
                            "batch_name": batch_name,
                            "input_bam_name": input_bam_name,
                            "positive_or_negative_loci": positive_or_negative_loci,
                            "locus": match.group(1),
                            "error message": match.group(0),
                        })

                # parse the log
                num_loci = int(re.search("Genotyping (\d+) loci", job_log["main"]).group(1))

                # example: "Elapsed (wall clock) time (h:mm:ss or m:ss): 7:53.62"
                wall_time_match = re.search("Elapsed [(]wall clock[)] time.*[)]: ([0-9:.]+)", job_log["main"])
                if not wall_time_match:
                    print("WARNING: couldn't find 'Elapsed (wall clock) time' in main log. Skipping...")
                    continue
                wall_time_string = wall_time_match.group(1)
                tool_run_time_seconds = sum(value * 60**i for i, value in enumerate(
                    [float(value) for value in wall_time_string.split(":")][::-1]  # parse time tokens in reverse order
                ))

                max_resident_set_match = re.search("Maximum resident set size [(]kbytes[)]: ([0-9]+)", job_log["main"])
                if not max_resident_set_match:
                    print("WARNING: couldn't find 'Maximum resident set size (kbytes)' line in main log. Skipping...")
                    continue
                max_resident_set_kbytes = int(max_resident_set_match.group(1))

                average_resident_set_match = re.search("Average resident set size [(]kbytes[)]: ([0-9]+)", job_log["main"])
                if not average_resident_set_match:
                    print("WARNING: couldn't find 'Average resident set size (kbytes)' line in main log. Skipping...")
                    continue
                average_resident_set_kbytes = float(average_resident_set_match.group(1))


                seconds_per_10k_loci = 10_000 * tool_run_time_seconds / num_loci
                minutes_per_10k_loci = seconds_per_10k_loci/60
                #print(f"num_loci: {num_loci:5d}   "
                #      f"total seconds: {tool_run_time_seconds:10.1f}   "
                #      f"minutes/10k loci: {minutes_per_10k_loci:10.1f} minutes")

                data_group_key = (tool_name, positive_or_negative_loci, coverage)
                output_row = {
                    "batch_name": batch_name,
                    "data_group_number": current_data_group_number[data_group_key],
                    "input_bam_name": input_bam_name,
                    "tool": tool_name,
                    "positive_or_negative_loci": positive_or_negative_loci,
                    "coverage": coverage,
                    "num_loci": num_loci,
                    "tool_run_time_seconds": tool_run_time_seconds,
                    "minutes_per_10k_loci": minutes_per_10k_loci,
                    "max_resident_set_kbytes": max_resident_set_kbytes,
                    "average_resident_set_kbytes": average_resident_set_kbytes,
                }
                output_rows.append(output_row)

                #print(" "*8, output_row)
                current_data_group_number[data_group_key] = (current_data_group_number[data_group_key] + 1) % number_of_data_groups


                job_stat = job_obj.status()
        except asyncio.TimeoutError as e:
            time.sleep(2)
            print("Retrying after timeout error...")
        else:
            break  # exit retry block
    else:
        print("All 3 retries failed... exiting.")
        break
bc.close()

#%%

output_dir = "./tool_comparison/hail_batch_pipelines"
for positive_or_negative in "positive", "negative":
    filename_prefix = "truth_set" if positive_or_negative == "positive" else "negative"
    output_filename = f"{output_dir}/{filename_prefix}_loci_that_cause_illumina_expansion_hunter_error.txt"
    counter = 0
    with open(output_filename, "wt") as f:
        for row in loci_with_expansion_hunter_errors:
            if row["positive_or_negative_loci"] == f"{positive_or_negative}_loci":
                counter += 1
                f.write(row["locus"] + "\n")
    print(f"Wrote {counter} rows to {output_filename}")


#%%

GROUP_BY_KEY = ['tool', 'positive_or_negative_loci', 'coverage', 'data_group_number']

df = pd.DataFrame(output_rows)
df2 = df.groupby(
    GROUP_BY_KEY
).aggregate({
    "num_loci": "sum",
    "tool_run_time_seconds": "sum",
    "max_resident_set_kbytes": "median",
    "average_resident_set_kbytes": "median",
}).reset_index()

df2.loc[:, "seconds_per_10k_loci"] = 10_000 * df2["tool_run_time_seconds"] / df2["num_loci"]
df2.loc[:, "minutes_per_10k_loci"] = df2["seconds_per_10k_loci"] / 60
df2.to_csv("./tool_comparison/STR_tool_timing.with_coverage.tsv", sep="\t", index=False, header=True)

#%%

#%%