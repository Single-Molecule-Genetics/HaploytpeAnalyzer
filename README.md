# Haplotype Analyzer

This tool is an addon to the Variant analyze. The tool extracts molecules wiht haplotypes(two or more variants on the same molecule), and classifies them through the same tier system used by the variant analyzer. This classification gives a confidence of what could be a true haplotype in a sequenced molecule versus sequenceing artefacts, and helps removing these artefacts which might appear as high quality variant, and pass the VarA filters. Removing these variants affect the allele frequencies which is recalculated through the tool and provided in a new Xlsx file.

## Dependencies
The tool works python 3.7 or higher,  and requires (pandas, openpyxl, xlsxwriter)

## Usage
A detailed description of all tools can be found on [Galaxy](http://usegalaxy.org), and on [JCU](https://invenio.nusl.cz/record/519820?ln=en) with all parameters, input and output files.

### Haplotype Analyser
Extracts molecules with haplotypes and outputs 3 Xlsx files, the haplotype analysis file, haplotypes tiers file, and an updated allele frequencies file.

**Input** 

**Dataset 1 (--SummaryFile):** XLSX summary file from the variant analyser output.

**Dataset 2 (--FreqFile):** XLSX variants frequencies file from the variant analyser.

**Output**

The outputs are, Xlsx file with extracted haplotypes Ref > Alt format, Xlsx file with extracted haplotypes tier number format, and Xlsx file with updated allele frequencies for the original VF file.

`$ python HaplotypeAnalysisV1.py -s $summaryfile -f $frequenciesfile --outputFile1 HaplotypeAnalysisV4.1.xlsx --outputFile2 Tier_AnalysisV4.1.xlsx --outputFile3 NewFreqV4.1.xlsx
