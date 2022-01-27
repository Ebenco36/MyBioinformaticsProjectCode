<h4>Basic Instructions</h4>

buildFiles(args.reference_genome_file)

python3 ischemicStroke.py --sra_run_file=SraRunTableIschemicStroke.txt --kallisto_quantify=yes

grep -e '^ENST[[:digit:]]\{11\}' ./gencode.v35.transcripts.fa > enst.txt

python3 ischemicStroke.py --sra_run_file=/Users/mac/Downloads/SraRunTableV2.txt

python3 ischemicStroke.py --sra_run_file=/Users/mac/Desktop/Bioinformatics/PythonScriptProject/SraRunTableBlood.txt --download_sra=yes

python3 ischemicStroke.py --sra_run_file=/Users/mac/Desktop/Bioinformatics/PythonScriptProject/SraRunTableBloodPheripheral.txt --download_sra=yes

python3 ischemicStroke.py --runAnalysis="yes" --result_path="sub_acute_control" --sra_run_file=SraRunTableIschemicStroke.txt --quantify_dir="quants/" --sample_split="1:5, 11:15" --split_group="SAS,NC", --sample_count="5, 5"
    
python3 ischemicStroke.py --runAnalysis="yes" --result_path="acute_control" --sra_run_file=SraRunTableIschemicStroke.txt --quantify_dir="quants/" --sample_split="6:10, 11:15" --split_group="AIS,NC", --sample_count="5, 5"

python3 ischemicStroke.py --runAnalysis="yes" --result_path="all_samples_with_control" --sra_run_file=SraRunTableIschemicStroke.txt --quantify_dir="quants/" --sample_split="1:5, 6:10, 11:15" --split_group="SAS, AIS, NC", --sample_count="5, 5, 5"
    
python3 ischemicStroke.py --runAnalysis="yes" --result_path="all_stroke_sample" --sra_run_file=SraRunTableIschemicStroke.txt --quantify_dir="quants/" --sample_split="1:5, 6:10" --split_group="SAS, AIS", --sample_count="5, 5"


python3 ischemicStroke.py --machineLearningAnalysis="yes" --input_file="../../../normalizedCounts.csv" --path_to_save="machineLearningReport" --component_n="3"  --split_group="SAS, AIS, NC", --sample_count="5, 5, 5" --graphFileName="result"

