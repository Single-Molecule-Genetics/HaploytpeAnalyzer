#!/usr/bin/env python

"""HaplotypeAnalyserV1.py
Author -- Shehab Moukbel Ali Aldawla
Contact -- shehab.moukbel@jku.at
Finds haplotypes within molecules and
classifies them based on quality driven from the variant analyser results.
=======  ==========  =================  ================================
Version  Date        Author             Description
4.1.0   2022-04-04  Shehab Moukbel    -
=======  ==========  =================  ================================
USAGE: python HaplotypeAnalyser.py --SummaryFile Variant analyser summary xlsx file --FreqFile Variant analyser
                                   variant's frequencies file --DcsBam DCS.bam
"""
from collections import Counter
import xlsxwriter
import argparse
import pandas as pd
import sys


def make_argparser():
    parser = argparse.ArgumentParser(description='Finds haplotypes withing  a library, and classifies them based'
                                                 'on tiers. It also updates the allele frequencies of the variants')

    parser.add_argument("-s", "--SummaryFile", type=str, required=True,
                        help='XLSX summary file from the variant analyser')
    parser.add_argument("-f", "--FreqFile", type=str, required=True,
                        help='XLSX variants frequencies file from the variant analyser')
    parser.add_argument('--outputFile1',
                        help='Output xlsx file with extracted haplotypes Ref > Alt format')
    parser.add_argument('--outputFile2',
                        help='Output xlsx file with extracted haplotypes tier number format')
    parser.add_argument('--outputFile3',
                        help='Output xlsx file with updated allele frequencies for the original VF file')

    return parser


def haplotype_analyser(argv):
    same_molecule = True
    parser = make_argparser()
    args = parser.parse_args(argv[1:])
    summary = args.SummaryFile
    frequency = args.FreqFile
    outfile1 = args.outputFile1
    outfile2 = args.outputFile2
    outfile3 = args.outputFile3

    # noinspection PyArgumentList
    df_s = pd.read_excel(summary)
    df_s['variant ID'] = df_s['variant ID'].fillna(',')
    df_s['in phase'] = df_s['in phase'].fillna(',')

    # noinspection PyArgumentList
    df_f = pd.read_excel(frequency, index_col=0)

    variants_freq_all = {var: df_f.loc[var]['AF (all tiers)'] for var in df_f.index.tolist()}
    variants_freq_gd = {var: df_f.loc[var]['AF (tiers 1.1-2.5)'] for var in df_f.index.tolist()}

    # Obtain the variants from the summary file and their tiers respectively depending on the tag they are on
    variants_tags_tier = {i: {} for i in df_s['variant ID']}
    for row in df_s.index.tolist():
        var = df_s.loc[row]['variant ID']
        tag = df_s.loc[row]['tag']
        tier = df_s.loc[row]['tier']
        variants_tags_tier[var][tag] = tier
    del variants_tags_tier[',']

    # get the linkages of each tag (variant id + in phase variants )
    # join in phase with variant ID and remove comas, white spaces
    linkages = [i.replace(' ', ',') for i in df_s['in phase']]
    linkages = [i + ',' + k for i, k in zip(linkages, df_s['variant ID'])]
    linkages = [list(set(i.split(','))) for i in linkages]
    df_s['linkages'] = linkages

    # get tags in a dictionary
    tags_haplotypes_all = {i: [] for i in df_s['tag']}
    for row in df_s.index.tolist():
        tag = df_s.loc[row]['tag']
        haplotype = df_s.loc[row]['linkages']  # add linkages to each tag
        tags_haplotypes_all[tag].append(haplotype)

    for key in list(tags_haplotypes_all):  # remove nan key
        if isinstance(key, float):
            del tags_haplotypes_all[key]

        # tags_haplotypes_all contains the tags with the haplotypes present on them and tags with only one variant
        # ( same molecule filter is not applied yet)
    tags_haplotypes_all = {tag: list(set(sum(tags_haplotypes_all[tag], []))) for tag in
                           tags_haplotypes_all}  # flatten the lists and get rid of duplicates
    tags_haplotypes_all = {tag: [var for var in tags_haplotypes_all[tag] if len(var) > 1] for tag in
                           tags_haplotypes_all}  # remove '' from the lists

    if same_molecule:  # Filter out those not on same molecule (sometimes issue that a variant has no tier)
        for tag in tags_haplotypes_all:
            try:
                haplotype = tags_haplotypes_all[tag]

                haplotype = [var for var in haplotype if tag in list(variants_tags_tier[var])]

                tags_haplotypes_all[tag] = haplotype
            except KeyError:
                haplotype = tags_haplotypes_all[tag]

                haplotype = [var for var in haplotype if
                             var in variants_tags_tier and tag in list(variants_tags_tier[var])]

                tags_haplotypes_all[tag] = haplotype

    swit_snps = [var for var in df_f.index.tolist() if
                 df_f.loc[var]['AF (tiers 1.1-2.5)'] >= 0.6]  # SNPS of 60% or higher to be removed from the haplotypes

    # All tiers sheet which will include all tags and their haplotypes regardless of their tiers
    all_tiers = {tag: tags_haplotypes_all[tag] for tag in tags_haplotypes_all if len(tags_haplotypes_all[tag]) > 1}

    min1_gd, min1_gd_f, final_hap, filtered_out, swit_haps = {}, {}, {}, {}, {}
    for tag in all_tiers:
        haplotype = all_tiers[tag]
        if any(var in swit_snps for var in haplotype):
            swit_haps[tag] = haplotype
            haplotype = [var for var in haplotype if var not in swit_snps]
            if len(haplotype) < 2:
                continue
        tiers = [variants_tags_tier[var][tag] for var in haplotype]

        if any(i < 3 for i in tiers):  # if any variant is from good tier then pass the haplotype to min1_gd
            min1_gd[tag] = haplotype
            filter_hap = [var for var in haplotype if
                          tiers[haplotype.index(var)] < 3]  # filter out variants of bad tiers  of the same haplotype

            if len(filter_hap) > 1:  # if the remaining haplotype is more than 1 variant after removing the bad tiers
                # then pass it to min1_gd_f
                min1_gd_f[tag] = filter_hap
                # final_hap
                counter = 0
                freqs = [variants_freq_gd[var] for var in filter_hap]
                for num in freqs:
                    if num < 0.01:
                        counter = counter + 1
                if counter < 2:
                    final_hap[tag] = filter_hap

    filtered_out = {tag: min1_gd_f[tag] for tag in min1_gd_f if tag not in final_hap}

    # V2 Updating frequency
    df_new = df_f[['cvrg (tiers 1.1-2.5)', 'AC alt (tiers 1.1-2.5)', 'AF (tiers 1.1-2.5)']]  # Obtain the original AF
    for var in df_new.index.tolist():  # Restrict the list to only those of high tier (forming haplotypes or lonely)
        if var in swit_snps or df_new.loc[var]['AF (tiers 1.1-2.5)'] == 0:
            df_new = df_new.drop(var)

    original_ac = {i: int(df_new.loc[i]['AC alt (tiers 1.1-2.5)']) for i in df_new.index.tolist() if
                   df_new.loc[i]['AC alt (tiers 1.1-2.5)'] != 0}  # put the variants in a dictionary with their AC
    # the occurrences of the variants in the filtered list ( variants that must have their AF/AC updated from the
    # original list)
    filtered_occ = sum(filtered_out.values(), [])
    filtered_occ = dict(Counter(filtered_occ))
    remained_ac = {}  # Results list (full variants)
    for i in original_ac:
        if i in filtered_occ:

            remained_ac[i] = original_ac[i] - filtered_occ[i]
        else:
            remained_ac[i] = original_ac[i]

    for variant in df_new.index.tolist():  # re updating the frequency dataframe and take only a sub part of relevant
        # columns in order to past it to a new Excel sheet
        if df_new.loc[variant]['AC alt (tiers 1.1-2.5)'] != remained_ac[variant] and remained_ac[variant] == 0:
            df_new = df_new.drop(variant)
        elif df_new.loc[variant]['AC alt (tiers 1.1-2.5)'] != remained_ac[variant] and remained_ac[variant] != 0:
            diff = (df_new.loc[variant]['AC alt (tiers 1.1-2.5)'] - remained_ac[variant])
            cvg = df_new.loc[variant]['cvrg (tiers 1.1-2.5)']
            df_new.at[variant, 'cvrg (tiers 1.1-2.5)'] = cvg - diff
            df_new.at[variant, 'AC alt (tiers 1.1-2.5)'] = remained_ac[variant]
            df_new.at[variant, 'AF (tiers 1.1-2.5)'] =\
                df_new.loc[variant]['AC alt (tiers 1.1-2.5)'] / df_new.loc[variant]['cvrg (tiers 1.1-2.5)']
    
    new_af_workbook = xlsxwriter.Workbook(outfile3)
    new_af_ws = new_af_workbook.add_worksheet("Allele frequencies")
    df_new_out =  df_new.reset_index()
    new_af_ws.write(0,0,'variant ID')
    new_af_ws.write(0,1,'cvrg (tiers 1.1-2.5)')
    new_af_ws.write(0,2,'AC alt (tiers 1.1-2.5)')
    new_af_ws.write(0,3,'AF (tiers 1.1-2.5)')
    for n_row_num, n_row_data in df_new_out.iterrows():
        for n_col_num, n_col_data in enumerate(row_data):
            new_af_ws.write(n_row_num+1, n_col_num, n_col_data)
    new_af_workbook.close()

    variants_freq_updated = {var: df_new.loc[var]['AF (tiers 1.1-2.5)'] for var in df_new.index.tolist()}
    tags_variants_tier = {tag: {var: variants_tags_tier[var][tag] for var in all_tiers[tag]} for tag in all_tiers}

    def analysis(analysis_sheet, analysis_writer, analysis_dic, analysis_freq_dic):
        try:
            ws = analysis_writer.add_worksheet(analysis_sheet)
            container = []

            worksheet_variants = [item for sublist in list(analysis_dic.values()) for item in sublist]
            worksheet_variants = sorted(list(set(worksheet_variants)))
            variants_freq = [analysis_freq_dic[i] for i in worksheet_variants]

            worksheet_columns = [i[:18] for i in worksheet_variants]

            for tag_analysis in list(analysis_dic):
                lst_temp = []

                for variant_a in worksheet_variants:
                    if variant_a in analysis_dic[tag_analysis]:
                        lst_temp.append(variant_a[-1])
                    else:
                        lst_temp.append('')
                container.append(lst_temp)

            df_result = pd.DataFrame(columns=worksheet_columns)
            df_result.loc['AF'] = ''

            for tag_analysis in list(analysis_dic):
                df_result.loc[tag_analysis] = container[(list(analysis_dic).index(tag_analysis))]
            
            df_result_temp = df_result.reset_index()
            for row_num, row_data in df_result_temp.iterrows():
                for col_num, col_data in enumerate(row_data):
                    ws.write(row_num + 1, col_num, col_data)


            red = analysis_writer.add_format({'bg_color': '#FF4F33'})
            green = analysis_writer.add_format({'bg_color': '#0FF235'})
            blue = analysis_writer.add_format({'bg_color': '#33A8FF'})
            orange = analysis_writer.add_format({'bg_color': '#FFB266'})
            pink = analysis_writer.add_format({'bg_color': '#FF99FF'})
            pink_dark = analysis_writer.add_format({'bg_color': '#FF00FF'})

            l_col = len(list(df_result))
            l_row = len(df_result.index.tolist())

            rotate_up = analysis_writer.add_format()
            rotate_up.set_rotation(90)

            rotate_angel = analysis_writer.add_format()
            rotate_angel.set_rotation(55)

            for one in range(len(variants_freq)):
                ws.write(0, one + 1, list(df_result)[one], rotate_angel)
                ws.write(1, one + 1, float(variants_freq[one]), rotate_up)

            ws.conditional_format(1, 1, l_row, l_col, {'type': 'text',
                                                                    'criteria': 'begins with',
                                                                    'value': 'A',
                                                                    'format': green})
            ws.conditional_format(1, 1, l_row, l_col, {'type': 'text',
                                                                    'criteria': 'begins with',
                                                                    'value': 'G',
                                                                    'format': orange})
            ws.conditional_format(1, 1, l_row, l_col, {'type': 'text',
                                                                    'criteria': 'begins with',
                                                                    'value': 'T',
                                                                    'format': red})
            ws.conditional_format(1, 1, l_row, l_col, {'type': 'text',
                                                                    'criteria': 'begins with',
                                                                    'value': 'C',
                                                                    'format': blue})
            ws.conditional_format(1, 1, 1, l_col, {'type': 'cell',
                                                                'criteria': 'between',
                                                                'minimum': 0.01,
                                                                'maximum': 0.4,
                                                                'format': pink})
            ws.conditional_format(1, 1, 1, l_col, {'type': 'cell',
                                                                'criteria': 'between',
                                                                'minimum': 0.41,
                                                                'maximum': 1,
                                                                'format': pink_dark})

            ws.set_column(0, 0, 30)
            ws.set_column(1, 9999, 2.33)
            ws.set_row(0, 95)
            ws.freeze_panes(2, 1)
        except ValueError:
            pass

    def analysis_tier(tier_sheet, tier_writer, tier_dic, variant_tier, t_dic_a):
        try:
            ws = tier_writer.add_worksheet(tier_sheet)
            tier_container = []
            tier_variants = [item for sublist in list(tier_dic.values()) for item in sublist]
            tier_variants = sorted(list(set(tier_variants)))
            tier_variants_freq = [t_dic_a[i] for i in tier_variants]
            tier_columns = [i[:18] for i in tier_variants]
           
            for tag_tier in list(tier_dic):
                tier_lst_temp = []

                for tier_variant in tier_variants:
                    if tier_variant in tier_dic[tag_tier]:
                        tier_lst_temp.append(str(variant_tier[tag_tier][tier_variant]))
                    else:
                        tier_lst_temp.append('')
                tier_container.append(tier_lst_temp)

            df_tier_result = pd.DataFrame(columns=tier_columns)
            df_tier_result.loc['AF'] = tier_variants_freq

            for tag_tier in list(tier_dic):
                df_tier_result.loc[tag_tier] = tier_container[(list(tier_dic).index(tag_tier))]


            df_tier_result_temp = df_tier_result.reset_index()

            for row_num, row_data in df_tier_result_temp.iterrows():
                for col_num, col_data in enumerate(row_data):
                    ws.write(row_num + 1, col_num, col_data)
            



            t_rotate_up = tier_writer.add_format()
            t_rotate_up.set_rotation(90)

            t_rotate_angel = tier_writer.add_format()
            t_rotate_angel.set_rotation(55)

            for one in range(len(tier_variants_freq)):
                ws.write(0, one + 1, list(df_tier_result)[one], t_rotate_angel)
                ws.write(1, one + 1, tier_variants_freq[one], t_rotate_up)

            l_col = len(list(df_tier_result))
            l_row = len(df_tier_result.index.tolist())
            t_red = tier_writer.add_format({'bg_color': '#FF4F33'})
            t_green = tier_writer.add_format({'bg_color': '#0FF235'})
            t_orange = tier_writer.add_format({'bg_color': '#FFB266'})
            t_pink = tier_writer.add_format({'bg_color': '#FF99FF'})
            t_pink_dark = tier_writer.add_format({'bg_color': '#FF00FF'})

            ws.conditional_format(2, 1, l_row, l_col, {'type': 'text',
                                                               'criteria': 'begins with',
                                                               'value': '1',
                                                               'format': t_green})
            ws.conditional_format(2, 1, l_row, l_col, {'type': 'text',
                                                               'criteria': 'begins with',
                                                               'value': '2',
                                                               'format': t_orange})
            ws.conditional_format(2, 1, l_row, l_col, {'type': 'text',
                                                               'criteria': 'begins with',
                                                               'value': '3',
                                                               'format': t_red})
            ws.conditional_format(2, 1, l_row, l_col, {'type': 'text',
                                                               'criteria': 'begins with',
                                                               'value': '4',
                                                               'format': t_red})
            ws.conditional_format(2, 1, l_row, l_col, {'type': 'text',
                                                               'criteria': 'begins with',
                                                               'value': '5',
                                                               'format': t_red})
            ws.conditional_format(2, 1, l_row, l_col, {'type': 'text',
                                                               'criteria': 'begins with',
                                                               'value': '6',
                                                               'format': t_red})
            ws.conditional_format(2, 1, l_row, l_col, {'type': 'text',
                                                               'criteria': 'begins with',
                                                               'value': '7',
                                                               'format': t_red})

            ws.conditional_format(1, 1, 1, l_col, {'type': 'cell',
                                                           'criteria': 'between',
                                                           'minimum': 0.01,
                                                           'maximum': 0.4,
                                                           'format': t_pink})
            ws.conditional_format(1, 1, 1, l_col, {'type': 'cell',
                                                           'criteria': 'between',
                                                           'minimum': 0.41,
                                                           'maximum': 1,
                                                           'format': t_pink_dark})

            ws.set_column(0, 0, 30)
            ws.set_column(1, 9999, 3.33)
            ws.set_row(0, 95)
            ws.freeze_panes(2, 1)
        except ValueError:
            pass

    sheet_names = ['All_tiers', 'Hap_1+', 'HapHQ', 'lowfreq 0r1', 'lowfreq >1', 'swit_haps']
    dic_lst = [all_tiers, min1_gd, min1_gd_f, final_hap, filtered_out, swit_haps]

    writer = xlsxwriter.Workbook(outfile1)
    writer_tier = xlsxwriter.Workbook(outfile2)


    for dictionary in dic_lst:  # to remove when a dictionary is empty

        if len(dictionary) < 1:
            sheet_names.pop(dic_lst.index(dictionary))
            dic_lst.pop(dic_lst.index(dictionary))

    for sheet, dic in zip(sheet_names, dic_lst):

        if sheet == 'lowfreq 0r1':
            f_dic_a = variants_freq_updated
        elif sheet == 'HapHQ' or sheet == 'lowfreq >1':
            f_dic_a = variants_freq_gd
        else:
            f_dic_a = variants_freq_all

        analysis(sheet, writer, dic, f_dic_a)
        analysis_tier(sheet, writer_tier, dic, tags_variants_tier, f_dic_a)

    writer.close()
    writer_tier.close()


if __name__ == '__main__':
    sys.exit(haplotype_analyser(sys.argv))
