import argparse
from numbers_utils import search, format_n, format_np
import os
import pandas as pd


def compute_num_failed_hg38_to_t2t_liftover(stepA_log_contents):
    # calculate # failed hg38 => T2T liftover
    print("\nLiftover failed for:")
    liftover_failed_IndelStraddlesMultipleIntevals = search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) IndelStraddlesMultipleIntevals", stepA_log_contents, type=int)

    liftover_failed_for_other_reason = 0
    liftover_failed_for_other_reason += search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) CannotLiftOver", stepA_log_contents, type=int)
    liftover_failed_for_other_reason += search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) MismatchedRefAllele", stepA_log_contents, type=int)
    liftover_failed_for_other_reason += search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) NoTarget", stepA_log_contents, type=int)

    total_STR_variants_before_validation_step = search(
        f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)


    print(f"{format_np(liftover_failed_IndelStraddlesMultipleIntevals, total_STR_variants_before_validation_step)} "
          f"STRs failed hg38 => T2T liftover due to IndelStraddlesMultipleIntevals error")
    print(f"{format_np(liftover_failed_for_other_reason, total_STR_variants_before_validation_step)} "
          f"STRs failed hg38 => T2T liftover for other reasons")

    return liftover_failed_for_other_reason


def compute_num_failed_comparison_with_t2t(stepA_log_contents):

    kept_het_variants = search(
        f"step4:STR[ ]+([0-9,]+) .*kept variants: heterozygous reference genotype", stepA_log_contents, type=int)

    kept_variants_insertions_that_match_adjacent_reference_sequence = search(
        f"step4:STR[ ]+([0-9,]+) .*kept variants: insertion matches the adjacent reference sequence", stepA_log_contents, type=int)

    total_variants_passed_t2t_comparison = search(
        f"step4:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

    assert total_variants_passed_t2t_comparison == kept_het_variants + kept_variants_insertions_that_match_adjacent_reference_sequence

    total_STR_variants_before_validation_step = search(
        f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

    print(f"{format_np(total_variants_passed_t2t_comparison, total_STR_variants_before_validation_step)} "
          "Total ")

    total_STR_variants_after_1st_liftover_step = search(
        f"step3:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

    num_failed = total_STR_variants_after_1st_liftover_step - total_variants_passed_t2t_comparison
    print(f"{format_np(num_failed, total_STR_variants_before_validation_step)} failed T2T comparison ")

    return num_failed


def compute_num_failed_t2t_to_hg38_liftover(stepA_log_contents):
    liftover_failed__no_target = search(f"step5:STR:failed-liftover:[ ]+([0-9,]+)[ ]+NoTarget", stepA_log_contents, type=int)
    total_STR_variants_before_validation_step = search(
        f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

    print(f"{format_np(liftover_failed__no_target, total_STR_variants_before_validation_step)} "
          "Failed T2T => hg38 liftover: No Target")

    return liftover_failed__no_target


def compute_num_had_different_position_after_hg38_to_t2t_to_hg38_liftover(stepA_log_contents):

    total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers = search(
        f"([0-9,]+) .*variants had a different position after hg38 => T2T => hg38", stepA_log_contents,
        expected_number_of_matches=1, type=int)
    print(f"{format_n(total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers)} "
          f"variants had a different position after hg38 => T2T => hg38 liftovers")

    return total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants-table", help="Path of variants table",  default="step2.STRs.variants.tsv.gz")
    parser.add_argument("--step-A-log",     help="Path of step A log file", default="step_A.log")
    parser.add_argument("--output-html",   help="Path of output table",    default="figures_and_tables/table2_validation.html")
    args = parser.parse_args()

    # check arg files exist
    for arg_name, arg_value in sorted(vars(args).items()):
        if arg_name.endswith("_table") or arg_name.endswith("_log"):
            if not os.path.exists(arg_value):
                parser.error(f"File does not exist: {arg_value}")

    # print arg values
    for arg_name, arg_value in sorted(vars(args).items()):
        print(f"{arg_name:30s} {arg_value}")


    with open(args.step_A_log, "rt") as f:
        stepA_log_contents = f.read()

    df_variants_before_validation = pd.read_table(args.variants_table)


    high_confidence_variants_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]variants", stepA_log_contents, type=int)
    high_confidence_alleles_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]alleles", stepA_log_contents, type=int)

    high_confidence_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_variants_in_syndip:,d}.*[)] SNV variants", stepA_log_contents).replace(",", ""))
    high_confidence_multiallelic_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_variants_in_syndip:,d}.*[)] multiallelic SNV variants", stepA_log_contents).replace(",", ""))
    high_confidence_indel_variants_in_syndip = high_confidence_variants_in_syndip - high_confidence_SNV_variants - high_confidence_multiallelic_SNV_variants

    high_confidence_INS_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] INS alleles", stepA_log_contents).replace(",", ""))
    high_confidence_DEL_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] DEL alleles", stepA_log_contents).replace(",", ""))
    high_confidence_indel_alleles_in_syndip = high_confidence_INS_alleles + high_confidence_DEL_alleles

    print("")
    print(f"{format_n(high_confidence_variants_in_syndip)} high-confidence variants in SynDip")
    print(f"{format_n(high_confidence_alleles_in_syndip)} high-confidence alleles in SynDip")

    print(f"{format_np(high_confidence_indel_variants_in_syndip, high_confidence_variants_in_syndip)} high-confidence INDEL variants in SynDip")
    print(f"{format_np(high_confidence_indel_alleles_in_syndip, high_confidence_alleles_in_syndip)} high-confidence INDEL alleles in SynDip")


    header = ["Step", "Description", "# failed", "% failed"]

    table_html = []
    table_html.append(f"<table>")
    table_html.append(f"<tr>{''.join(f'<th nowrap>{h}</th>' for h in header)}</tr>")

    # row 0
    table_html.append("<tr>")
    table_html.append("""<td><i></i></td>""")
    table_html.append(f"""<td>Start with all {len(df_variants_before_validation):,d} TR variants in SynDip high-confidence regions</td>""")
    table_html.append(f"""<td></td>""")
    table_html.append(f"""<td></td>""")
    table_html.append("</tr>")

    # row 1
    num_variants_failed_hg38_to_t2t_liftover = compute_num_failed_hg38_to_t2t_liftover(stepA_log_contents)
    table_html.append("<tr>")
    table_html.append("""<td><i>1</i></td>""")
    table_html.append("""<td>Liftover TR variants from hg38 ⇒ T2T</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_failed_hg38_to_t2t_liftover:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_failed_hg38_to_t2t_liftover/len(df_variants_before_validation):.1%}</td>""")
    table_html.append("</tr>")

    # row 2
    num_variants_failed_comparison_with_t2t = compute_num_failed_comparison_with_t2t(stepA_log_contents)
    table_html.append("<tr>")
    table_html.append("""<td><i>2</i></td>""")
    table_html.append("""<td>Check that at least one allele matches the T2T reference sequence</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_failed_comparison_with_t2t:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_failed_comparison_with_t2t/len(df_variants_before_validation):.1%}</td>""")
    table_html.append("</tr>")

    # row 3
    num_variants_failed_t2t_to_hg38_liftover = compute_num_failed_t2t_to_hg38_liftover(stepA_log_contents)
    table_html.append("<tr>")
    table_html.append("""<td><i>3</i></td>""")
    table_html.append("""<td>Lift variants back over from T2T ⇒ hg38</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_failed_t2t_to_hg38_liftover:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_failed_t2t_to_hg38_liftover/len(df_variants_before_validation):.1%}</td>""")
    table_html.append("</tr>")

    # row 4
    num_variants_had_different_position_after_hg38_to_t2t_to_hg38_liftover = \
        compute_num_had_different_position_after_hg38_to_t2t_to_hg38_liftover(stepA_log_contents)
    table_html.append("<tr>")
    table_html.append("""<td><i>4</i></td>""")
    table_html.append("""<td>Check if the variant's position changed after hg38 ⇒ T2T ⇒ hg38 liftover</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_had_different_position_after_hg38_to_t2t_to_hg38_liftover:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{num_variants_had_different_position_after_hg38_to_t2t_to_hg38_liftover/len(df_variants_before_validation):.1%}</td>""")
    table_html.append("</tr>")

    table_html.append("</table>")

    with open("figures_and_tables/table_template.html", "rt") as f, open(args.output_html, "wt") as fo:
        table_html_template = f.read()
        table_html_template = table_html_template % {
            "title": "Validation",
            "table": "\n".join(table_html),
        }
        fo.write(table_html_template)


if __name__ == "__main__":
    main()

