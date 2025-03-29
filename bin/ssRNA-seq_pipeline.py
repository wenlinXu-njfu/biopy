#!/usr/bin/env python
"""
File: ssRNA-seq_pipeline.py
Description: Strand specific (dUTP constructed library) pair end RNA-seq analysis pipeline (including lncRNA and its target gene prediction).
Date: 2023/2/21
Author: xuwenlin
E-mail: wenlinxu.njfu@outlook.com
"""
from io import TextIOWrapper
from collections import defaultdict
from os import makedirs, system
from os.path import abspath
from shutil import which
from yaml import safe_load
from schema import Schema, And, Use, SchemaError
import click
from pybioinformatic import (
    parse_sample_info,
    MergeSamples,
    StrandSpecificRNASeqAnalyser,
    LncRNAPredictor,
    LncRNAClassification,
    TaskManager,
    Displayer
)
displayer = Displayer(__file__.split('/')[-1], version='0.6.0')


def check_dependency():
    software_list = ['fastp', 'hisat2', 'gffread', 'stringtie', 'cuffcompare', 'featureCounts',
                     'CNCI.py', 'CPC2.py', 'PLEK', 'pfam_scan.pl', 'bedtools', 'R']
    flag = []
    click.echo('\033[36mCheck dependency.\033[0m', err=True)
    for software in software_list:
        path = which(software)
        if path:
            click.echo(f'{software}: {path}', err=True)
            flag.append(True)
        else:
            click.echo(f'{software}: command not found', err=True)
            flag.append(False)
    tkm = TaskManager(num_processing=1)
    cmd1 = '''R CMD Rscript -e 'if(requireNamespace("DESeq2", quietly = TRUE)) {print("True")} else {print("False")}' '''
    stdout1 = tkm.echo_and_exec_cmd(cmd=cmd1, show_cmd=False)
    if 'True' in stdout1:
        flag.append(True)
    else:
        flag.append(False)
        click.echo('\033[31mR package DESeq2 has not been installed.\033[0m', err=True)
    cmd2 = '''R CMD Rscript -e 'if(requireNamespace("clusterProfiler", quietly = TRUE)) {print("True")} else {print("False")}' '''
    stdout2 = tkm.echo_and_exec_cmd(cmd=cmd2, show_cmd=False)
    if 'True' in stdout2:
        flag.append(True)
    else:
        flag.append(False)
        click.echo('\033[31mR package clusterProfiler has not been installed.\033[0m', err=True)
    if not all(flag):
        exit()


def check_config(yaml_file: TextIOWrapper):
    click.echo('\033[36mCheck params.\033[0m', err=True)
    config_schema = Schema({
        "input": {
            "sample_info": str,
            "ref_genome": str,
            "ref_genome_gff": str,
            "enrich_anno_file": str
        },
        "output": {"dir": str},
        "global_params": {
            "num_threads": And(Use(int), lambda x: 0 < x, error="num_threads must be positive integer."),
            "num_processing": And(Use(int), lambda x: 0 < x, error="num_processing must be positive integer.")
        },
        "featureCounts_params": {
            "feature_type": click.Choice(['gene', 'mRNA', 'transcript', 'five_prime_UTR' 'CDS', 'three_prime_UTR', 'exon']),
            "mate_feature": click.Choice(['ID', 'Name', 'gene_id', 'transcript_id'])
        },
        "CNCI_params": {
            "CNCI_module": str
        },
        "pfamscan_params": {
            "pfamscan_database": str
        },
        "lncRNA_target_prediction_params": {
            "lncRNA_min_exp": float,
            "mRNA_min_exp": float,
            "r": And(Use(float), lambda x: 0 <= x <= 1, error="r must between 0 and 1."),
            "FDR": And(Use(float), lambda x: 0 <= x <= 1, error="FDR must between 0 and 1."),
            "q_value": And(Use(float), lambda x: 0 <= x <= 1, error="q_value must between 0 and 1."),
            "distance": And(Use(int), lambda x: 0 < x, error="distance must be positive integer.")
        },
        "DESeq2_params": {
            "padj": And(Use(float), lambda x: 0 <= x <= 1, error="padj must between 0 and 1."),
            "log2FoldChange": And(Use(float), lambda x: 0 < x, error="log2FoldChange must be greater than 0.")
        },
        "clusterProfiler_params": {
            "pvalueCutoff": And(Use(float), lambda x: 0 <= x <= 1, error="pvalueCutoff must between 0 and 1."),
            "pAdjustMethod": click.Choice(["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]),
            "qvalueCutoff": And(Use(float), lambda x: 0 <= x <= 1, error="qvalueCutoff must between 0 and 1.")
        }
    })

    config = safe_load(yaml_file)
    try:
        config_schema.validate(config)
        click.echo("\033[32mConfiguration verification passed.\n\033[0m", err=True)
    except SchemaError as e:
        click.echo(f"\033[31mConfig error: {e}\033[0m", err=True)
        exit()
    else:
        return config


def main(config: TextIOWrapper):
    check_dependency()
    config = check_config(config)

    # input
    sample_info = config['input']['sample_info']
    genome = config['input']['ref_genome']
    gff = config['input']['ref_genome_gff']
    kegg_anno_file = config['input']['enrich_anno_file']

    # output
    output_path = config['output']['dir']

    # global params
    num_threads = config['global_params']['num_threads']
    num_processing = config['global_params']['num_processing']

    # featureCounts params
    feature_type = config['featureCounts_params']['feature_type']
    count_unit = config['featureCounts_params']['mate_feature']

    # CNCI params
    module = config['CNCI_params']['CNCI_module']

    # pfamscan params
    pfamscan_database = config['pfamscan_params']['pfamscan_database']

    # lncRNA target prediction params
    lncRNA_min_exp = config['lncRNA_target_prediction_params']['lncRNA_min_exp']
    mRNA_min_exp = config['lncRNA_target_prediction_params']['mRNA_min_exp']
    r = config['lncRNA_target_prediction_params']['r']
    FDR = config['lncRNA_target_prediction_params']['FDR']
    q_value = config['lncRNA_target_prediction_params']['q_value']
    distance = config['lncRNA_target_prediction_params']['distance']

    # DESeq2 params
    padj = config['DESeq2_params']['padj']
    log2FoldChange = config['DESeq2_params']['log2FoldChange']

    # clusterProfiler params
    pvalueCutoff = config['clusterProfiler_params']['pvalueCutoff']
    pAdjustMethod = config['clusterProfiler_params']['pAdjustMethod']
    qvalueCutoff = config['clusterProfiler_params']['qvalueCutoff']

    exec_cmds = which('exec_cmds')
    file_split = which('file_split')
    featureCounts = which('featureCounts')
    featureCounts_helper = which('featureCounts_helper')
    joint = which('joint')
    ORF_finder = which('ORF_finder')
    plot = which('plot')
    python = which('python')
    Rscript = which('Rscript')
    seqkit = which('seqkit')

    output_path = abspath(output_path)
    storer = MergeSamples(set(), set())
    sample_info_dict = parse_sample_info(sample_info=sample_info)
    script_template_path = '/'.join(__file__.split('/')[:-2]) + '/script_template'
    comparative_combination = defaultdict(lambda: defaultdict(list))  # {C1: {C: [sample1, sample2, ...], T: [sample3, sample4, ...}, ...}

    # single sample analyse
    for sample_name, l in sample_info_dict.items():
        cc_list = l[2].split(';')
        for i in cc_list:
            cc_name, sample_class = i.split(':')
            comparative_combination[cc_name][sample_class].append(f'{sample_name}.ht2.sort.bam')
        makedirs(f'{output_path}/shell/normal', exist_ok=True)
        analyser = StrandSpecificRNASeqAnalyser(
            read1=l[0],
            read2=l[1],
            ref_gff=gff,
            ref_genome=genome,
            output_path=output_path,
            num_threads=num_threads,
            sample_name=sample_name
        )
        with open(f'{output_path}/shell/normal/{sample_name}.sh', 'w') as f:
            fastp = analyser.run_fastp(q_value=20, fastp_exe='fastp')
            hisat2 = analyser.run_hisat2(hisat2_exe='hisat2', samtools_exe='samtools', storer=storer)
            stringtie = analyser.run_stringtie(stringtie_exec='stringtie', storer=storer)
            cmds = [fastp, hisat2, stringtie]
            f.write('\n'.join(cmds))
    system(
        f'for i in `ls {output_path}/shell/normal`;'
        f'do echo sh {output_path}/shell/normal/$i;'
        f'done > {output_path}/shell/merge_normal.sh'
    )

    # merge results of each sample
    featureCounts_num_threads = num_threads * num_processing if num_threads * num_processing <= 64 else 64
    a = StrandSpecificRNASeqAnalyser(
        read1='',
        read2='',
        ref_genome=genome,
        output_path=output_path,
        num_threads=featureCounts_num_threads,
        ref_gff=gff
    )

    stringtie_merge = a.run_stringtie_merge(
        gtf_list=list(storer.stringtie_gtf),
        m=200,
        c=1.0,
        F=0.5,
        g=500,
        stringtie_exec='stringtie'
    )

    cuffcompare = a.run_cuffcompare(cuffcompare_exe='cuffcompare', gffread_exe='gffread')

    anno_file = 'gtf' if feature_type == 'exon' else 'gff'
    all_featureCounts = a.run_featureCounts(
        anno_file=anno_file,
        bam_list=list(storer.hisat2_bam),
        feature_type=feature_type,
        count_unit=count_unit,
        featureCounts_exec='featureCounts'
    )

    # lncRNA prediction
    makedirs(f'{output_path}/05.lncRNA_prediction', exist_ok=True)
    lp = LncRNAPredictor(
        nucl_fasta_file=f'{output_path}/03.assembly/novel_transcript.fa',
        output_path=f'{output_path}/05.lncRNA_prediction',
        pep_fasta_file=None,
        num_thread=num_threads,
        module=module
    )
    CNCI_cmd = lp.run_CNCI(CNCI_exec='CNCI.py')
    CPC2_cmd = lp.run_CPC2(CPC2_exec='CPC2.py')
    PLEK_cmd = lp.run_PLEK(PLEK_exec='PLEK')
    pfamscan_helper = which('pfamscan_helper')
    pfamscan_cmd = (
        f'{pfamscan_helper} batch '
        f'-i {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_input '
        f'-d {pfamscan_database} '
        f'-o {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_out '
        f'-n {num_threads * num_processing}'
    )
    lncRNA_prediction_cmd = f'{CNCI_cmd}\n{CPC2_cmd}\n{PLEK_cmd}\n{pfamscan_cmd}'
    lncRNA_results = lp.merge_results(
        CNCI_results=f'{output_path}/05.lncRNA_prediction/CNCI/CNCI.index',
        CPC2_results=f'{output_path}/05.lncRNA_prediction/CPC2/CPC2.txt',
        PLEK_results=f'{output_path}/05.lncRNA_prediction/PLEK/PLEK.xls',
        PfamScan_results=f'{output_path}/05.lncRNA_prediction/PfamScan/pfamscan_out/all_results.xls',
        seqkit_exec=seqkit
    )
    with open(f'{output_path}/shell/lncRNA_prediction.sh', 'w') as o:
        o.write(lncRNA_prediction_cmd)

    # lncRNA target gene prediction
    makedirs(f'{output_path}/06.lncRNA_target_prediction', exist_ok=True)
    with open(f'{script_template_path}/lncRNA_target_prediction') as f, \
        open(f'{output_path}/shell/lncRNA_target_prediction.py', 'w') as o:
        target_prediction_script = f.read().format(
            out_path=output_path,
            lncRNA_min_exp=lncRNA_min_exp,
            mRNA_min_exp=mRNA_min_exp,
            r=r,
            FDR=FDR,
            q_value=q_value,
            distance=distance,
            num_processing=num_processing * num_threads
        )
        o.write(target_prediction_script)

    # lncRNA classification
    lc = LncRNAClassification(
        mRNA_gff_file=gff,
        lncRNA_gtf_file=f'{output_path}/05.lncRNA_prediction/lncRNA.gtf',
        out_dir=f'{output_path}/07.lncRNA_classification'
    )
    classify_script = lc.classification()
    sense = (r'''cut -f 4 %s | cut -d' ' -f 4 | sed 's/"//g;s/;//' | awk '{print $0"\tsense"}' ''' %
             f'{output_path}/07.lncRNA_classification/sense.bed')
    antisense = (r'''cut -f 4 %s | cut -d' ' -f 4 | sed 's/"//g;s/;//' | awk '{print $0"\tantisense"}' ''' %
                 f'{output_path}/07.lncRNA_classification/antisense.bed')
    intronic = (r'''cut -f 4 %s | cut -d' ' -f 4 | sed 's/"//g;s/;//' | awk '{print $0"\tintronic"}' ''' %
                f'{output_path}/07.lncRNA_classification/intronic.bed')
    intergenic = (r'''cut -f 4 %s | cut -d' ' -f 4 | sed 's/"//g;s/;//' | awk '{print $0"\tintergenic"}' ''' %
                  f'{output_path}/07.lncRNA_classification/intergenic.bed')
    cmd = (fr"cat <({sense}) <({antisense}) <({intronic}) <({intergenic}) | sort -uV | sed '1iLncRNA_id\tLncRNA_type' | "
           fr"{joint} -i {output_path}/06.lncRNA_target_prediction/co_loc.xls "
           fr"-I - -o {output_path}/07.lncRNA_classification/co_loc.xls")
    with open(f'{output_path}/shell/lncRNA_classification.sh', 'w') as o:
        click.echo(classify_script, o)
        click.echo(cmd, o)

    # mRNA differential expression and enrichment analyse
    makedirs(f'{output_path}/shell/DESeq2/mRNA', exist_ok=True)
    makedirs(f'{output_path}/shell/clusterProfiler/mRNA', exist_ok=True)
    makedirs(f'{output_path}/shell/DE_enrich', exist_ok=True)
    with open(f'{script_template_path}/DESeq2') as f1, open(f'{script_template_path}/clusterProfiler') as f2:
        raw_DESeq2_script = f1.read()
        raw_clusterProfiler_scrpit = f2.read()
        for cc_name, d in comparative_combination.items():
            makedirs(f'{output_path}/08.mRNA_differential_expression_analysis/{cc_name}', exist_ok=True)
            makedirs(f'{output_path}/10.DEmRNA_enrich/{cc_name}', exist_ok=True)
            with open(f'{output_path}/shell/DE_enrich/DEmRNA_{cc_name}.sh', 'w') as o1:
                control_samples = ','.join(d['C'])
                num_control_sample = len(d['C'])
                treat_samples = ','.join(d['T'])
                num_treat_sample = len(d['T'])
                cmd = r'''awk -F'\t' -v cols="Geneid,{control_samples},{treat_samples}" '''.format(
                    control_samples=control_samples,
                    treat_samples=treat_samples
                )
                cmd += r''''BEGIN { split(cols, colarr, /,/); } NR==1 { for(i=1;i<=NF;i++) colname[$i]=i; n=0;for(j=1;j<=length(colarr);j++) if (colarr[j] in colname) indices[++n]=colname[colarr[j]]; } { for(k=1;k<=n;k++) printf "%s%s", $indices[k], (k==n?"\n":FS) }' '''
                cmd += '''{reads_count_path} > {out_file}\n'''.format(
                    reads_count_path=f'{output_path}/06.lncRNA_target_prediction/target_exp/reads.count.fc.xls',
                    out_file=f'{output_path}/08.mRNA_differential_expression_analysis/{cc_name}/reads.count.fc.xls'
                )
                cmd += (
                    f'{Rscript} {output_path}/shell/DESeq2/mRNA/{cc_name}.R > {output_path}/08.mRNA_differential_expression_analysis/{cc_name}/DESeq2.log 2>&1\n'
                    f'{plot} volcano -i {output_path}/08.mRNA_differential_expression_analysis/{cc_name}/volcano_data.txt -l 1.5 -p 0.05 -o {output_path}/08.mRNA_differential_expression_analysis/{cc_name}/volcano.pdf\n\n'
                )
                cmd += (
                    f"sed '1d' {output_path}/08.mRNA_differential_expression_analysis/{cc_name}/DESeq2.csv | cut -d',' -f 1 > {output_path}/10.DEmRNA_enrich/{cc_name}/DEmRNA.lst\n"
                    f"{Rscript} {output_path}/shell/clusterProfiler/mRNA/{cc_name}.R > {output_path}/10.DEmRNA_enrich/{cc_name}/clusterProfiler.log 2>&1"
                )
                o1.write(cmd)
                with open(f'{output_path}/shell/DESeq2/mRNA/{cc_name}.R', 'w') as o2:
                    DESeq2_script = raw_DESeq2_script.format(
                        reads_count_file=f'{output_path}/08.mRNA_differential_expression_analysis/{cc_name}/reads.count.fc.xls',
                        num_control_repeat=num_control_sample,
                        num_treat_repeat=num_treat_sample,
                        out_path=f'{output_path}/08.mRNA_differential_expression_analysis/{cc_name}',
                        padj=padj,
                        log2FoldChange=log2FoldChange
                    )
                    o2.write(DESeq2_script)
                with open(f'{output_path}/shell/clusterProfiler/mRNA/{cc_name}.R', 'w') as o3:
                    clusterProfiler_script = raw_clusterProfiler_scrpit.format(
                        out_path=f'{output_path}/10.DEmRNA_enrich/{cc_name}',
                        anno_file=kegg_anno_file,
                        DEgenes=f'{output_path}/10.DEmRNA_enrich/{cc_name}/DEmRNA.lst',
                        pvalueCutoff=pvalueCutoff,
                        pAdjustMethod=pAdjustMethod,
                        qvalueCutoff=qvalueCutoff
                    )
                    o3.write(clusterProfiler_script)

    # lncRNA differential expression and target enrichment analyse
    makedirs(f'{output_path}/shell/DESeq2/lncRNA', exist_ok=True)
    makedirs(f'{output_path}/shell/clusterProfiler/co_loc', exist_ok=True)
    makedirs(f'{output_path}/shell/clusterProfiler/co_exp', exist_ok=True)
    for cc_name, d in comparative_combination.items():
        makedirs(f'{output_path}/09.lncRNA_differential_expression_analysis/{cc_name}', exist_ok=True)
        makedirs(f'{output_path}/11.DElncRNA_target_enrich/{cc_name}/co_loc', exist_ok=True)
        makedirs(f'{output_path}/11.DElncRNA_target_enrich/{cc_name}/co_exp', exist_ok=True)
        with open(f'{output_path}/shell/DE_enrich/DElncRNA_{cc_name}.sh', 'w') as o1:
            control_samples = ','.join(d['C'])
            num_control_sample = len(d['C'])
            treat_samples = ','.join(d['T'])
            num_treat_sample = len(d['T'])
            cmd = r'''awk -F'\t' -v cols="Geneid,{control_samples},{treat_samples}" '''.format(
                control_samples=control_samples,
                treat_samples=treat_samples
            )
            cmd += r''''BEGIN { split(cols, colarr, /,/); } NR==1 { for(i=1;i<=NF;i++) colname[$i]=i; n=0;for(j=1;j<=length(colarr);j++) if (colarr[j] in colname) indices[++n]=colname[colarr[j]]; } { for(k=1;k<=n;k++) printf "%s%s", $indices[k], (k==n?"\n":FS) }' '''
            cmd += '''{reads_count_path} > {out_file}\n'''.format(
                reads_count_path=f'{output_path}/06.lncRNA_target_prediction/lncRNA_exp/reads.count.fc.xls',
                out_file=f'{output_path}/09.lncRNA_differential_expression_analysis/{cc_name}/reads.count.fc.xls'
            )
            cmd += (
                f'{Rscript} {output_path}/shell/DESeq2/lncRNA/{cc_name}.R > {output_path}/09.lncRNA_differential_expression_analysis/{cc_name}/DESeq2.log 2>&1\n'
                f'{plot} volcano -i {output_path}/09.lncRNA_differential_expression_analysis/{cc_name}/volcano_data.txt -l 1.5 -p 0.05 -o {output_path}/09.lncRNA_differential_expression_analysis/{cc_name}/volcano.pdf\n\n'
            )
            cmd += (
                f"sed '1d' {output_path}/09.lncRNA_differential_expression_analysis/{cc_name}/DESeq2.csv | cut -d',' -f 1 | grep -Fw -f - {output_path}/06.lncRNA_target_prediction/co_loc.xls | cut -f 2 | sort -uV > {output_path}/11.DElncRNA_target_enrich/{cc_name}/co_loc/target.lst\n"
                f"{Rscript} {output_path}/shell/clusterProfiler/co_loc/{cc_name}.R > {output_path}/11.DElncRNA_target_enrich/{cc_name}/co_loc/clusterProfiler.log 2>&1\n\n"
            )
            cmd += (
                f"sed '1d' {output_path}/09.lncRNA_differential_expression_analysis/{cc_name}/DESeq2.csv | cut -d',' -f 1 | grep -Fw -f - {output_path}/06.lncRNA_target_prediction/filter_co_exp.xls | cut -f 2 | sort -uV > {output_path}/11.DElncRNA_target_enrich/{cc_name}/co_exp/target.lst\n"
                f"{Rscript} {output_path}/shell/clusterProfiler/co_exp/{cc_name}.R > {output_path}/11.DElncRNA_target_enrich/{cc_name}/co_exp/clusterProfiler.log 2>&1\n\n"
            )
            o1.write(cmd)
            with open(f'{output_path}/shell/DESeq2/lncRNA/{cc_name}.R', 'w') as o2:
                DESeq2_script = raw_DESeq2_script.format(
                    reads_count_file=f'{output_path}/09.lncRNA_differential_expression_analysis/{cc_name}/reads.count.fc.xls',
                    num_control_repeat=num_control_sample,
                    num_treat_repeat=num_treat_sample,
                    out_path=f'{output_path}/09.lncRNA_differential_expression_analysis/{cc_name}',
                    padj=padj,
                    log2FoldChange=log2FoldChange
                )
                o2.write(DESeq2_script)
            with open(f'{output_path}/shell/clusterProfiler/co_loc/{cc_name}.R', 'w') as o3:
                clusterProfiler_script = raw_clusterProfiler_scrpit.format(
                    out_path=f'{output_path}/11.DElncRNA_target_enrich/{cc_name}/co_loc',
                    anno_file=kegg_anno_file,
                    DEgenes=f'{output_path}/11.DElncRNA_target_enrich/{cc_name}/co_loc/target.lst',
                    pvalueCutoff=pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    qvalueCutoff=qvalueCutoff
                )
                o3.write(clusterProfiler_script)
            with open(f'{output_path}/shell/clusterProfiler/co_exp/{cc_name}.R', 'w') as o4:
                clusterProfiler_script = raw_clusterProfiler_scrpit.format(
                    out_path=f'{output_path}/11.DElncRNA_target_enrich/{cc_name}/co_exp',
                    anno_file=kegg_anno_file,
                    DEgenes=f'{output_path}/11.DElncRNA_target_enrich/{cc_name}/co_exp/target.lst',
                    pvalueCutoff=pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    qvalueCutoff=qvalueCutoff
                )
                o4.write(clusterProfiler_script)

    # write all step commands
    with open(f'{output_path}/shell/All_step.sh', 'w') as o:
        mapping_rate_summary = r'''grep over %s/02.mapping/*/*.log | awk -F'/' '{print $NF}' | cut -d' ' -f 1 | sed 's/:/\t/' | sed 's/.ht2.log//' | sort -uV > %s/02.mapping/mapping_rate_summary.xls''' % (output_path, output_path)
        CNCI_results = r'''awk '{if($3 =="noncoding"){print $1}}' %s/05.lncRNA_prediction/CNCI/CNCI.index | sort -uV''' % output_path
        CPC2_results = r'''awk '{if($9 == "noncoding"){print $1}}' %s/05.lncRNA_prediction/CPC2/CPC2.txt | sort -uV''' % output_path
        PLEK_results = r'''awk '{if($1 == "Non-coding"){print $3}}' %s/05.lncRNA_prediction/PLEK/PLEK.xls | sed 's/>//' | sort -uV''' % output_path
        PfamScan_results = r'''cut -f1 %s/05.lncRNA_prediction/PfamScan/pfamscan_out/all_results.xls | grep -vFw -f - <(grep '>' %s/03.assembly/novel_transcript.fa | sed 's/>//') | cut -d' ' -f 1 | sort -uV''' % (output_path, output_path)
        plot_venn = f'''{plot} venn -g CNCI,CPC2,PfamScan,PLEK <({CNCI_results}) <({CPC2_results}) <({PfamScan_results}) <({PLEK_results}) -o {output_path}/05.lncRNA_prediction/venn.pdf'''
        create_lncRNA_gtf = (
            f'grep "^>" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -Fwf - {output_path}/03.assembly/All.gtf > {output_path}/05.lncRNA_prediction/lncRNA.gtf'
        )
        awk1 = r'''awk '{if(match($0, /transcript_id "[a-zA-Z0-9]*.*[a-zA-Z0-9]*"/)) print substr($0, RSTART, RLENGTH)}' %s''' % f'{output_path}/05.lncRNA_prediction/lncRNA.gtf'
        awk2 = r'''awk '{print $2"\t"$1}' '''
        lncRNA_exon_count = f'''{awk1} | cut -d' ' -f 2 | sed 's/"//g' | sort | uniq -c | {awk2}> {output_path}/05.lncRNA_prediction/lncRNA_exon_count.xls'''
        makedirs(f'{output_path}/06.lncRNA_target_prediction/lncRNA_exp', exist_ok=True)
        makedirs(f'{output_path}/06.lncRNA_target_prediction/target_exp', exist_ok=True)
        lncRNA_exp = (
            f'{featureCounts} -t {feature_type} -g {count_unit} -fMO -p --countReadPairs -s 2 -T {featureCounts_num_threads} '
            f'-a {output_path}/05.lncRNA_prediction/lncRNA.gtf '
            f'-o {output_path}/06.lncRNA_target_prediction/lncRNA_exp/lncRNA.fc.xls '
            f'{output_path}/02.mapping/*/*.ht2.sort.bam '
            f'2> {output_path}/06.lncRNA_target_prediction/lncRNA_exp/lncRNA.fc.log\n'
            f'{featureCounts_helper} normalization '
            f'-i {output_path}/06.lncRNA_target_prediction/lncRNA_exp/lncRNA.fc.xls '
            f'-o {output_path}/06.lncRNA_target_prediction/lncRNA_exp'
        )
        create_target_gtf = (
            f'grep "^>" {output_path}/05.lncRNA_prediction/lncRNA.fa | '
            f'sed "s/>//;s/ .*//" | '
            f'grep -vFwf - {output_path}/03.assembly/All.gtf > {output_path}/06.lncRNA_target_prediction/target.gtf'
        )
        target_exp = (
            f'{featureCounts} -t {feature_type} -g {count_unit} -fMO -p --countReadPairs -s 2 -T {featureCounts_num_threads} '
            f'-a {output_path}/06.lncRNA_target_prediction/target.gtf '
            f'-o {output_path}/06.lncRNA_target_prediction/target_exp/target.fc.xls '
            f'{output_path}/02.mapping/*/*.ht2.sort.bam '
            f'2> {output_path}/06.lncRNA_target_prediction/target_exp/target.fc.log\n'
            f'{featureCounts_helper} normalization '
            f'-i {output_path}/06.lncRNA_target_prediction/target_exp/target.fc.xls '
            f'-o {output_path}/06.lncRNA_target_prediction/target_exp'
        )
        joint = (fr"{joint} -i {output_path}/06.lncRNA_target_prediction/co_loc.xls -lc 6 "
                 fr"-I <(sed '1iLncRNA_id\tLncRNA_exon_count' {output_path}/05.lncRNA_prediction/lncRNA_exon_count.xls) "
                 fr"-o {output_path}/06.lncRNA_target_prediction/co_loc.xls")
        run_DE_enrich = (
            f'for i in `ls {output_path}/shell/DE_enrich`; do echo sh {output_path}/shell/DE_enrich/$i; done | {exec_cmds} -f - -n {num_threads}'
        )
        cmd = (
            f'#!/bin/bash\n\n'
            f'{exec_cmds} -f {output_path}/shell/merge_normal.sh -n {num_processing}\n\n'
            f'{mapping_rate_summary}\n\n'
            f'{stringtie_merge}\n\n'
            f'{cuffcompare}\n\n'
            f'{all_featureCounts}\n\n'
            f'{featureCounts_helper} normalization '
            f'-i {output_path}/04.expression/featureCounts.xls '
            f'-o {output_path}/04.expression\n\n'
            f'{ORF_finder} -l 30 -F -log {output_path}/03.assembly/ORF.log '
            f'-n {num_threads * num_processing} '
            f'-o {output_path}/03.assembly {output_path}/03.assembly/novel_transcript.fa\n\n'
            f'{file_split} -n 200 '
            f'-i {output_path}/03.assembly/novel_transcript_pep.fa '
            f'-o {output_path}/05.lncRNA_prediction/PfamScan/pfamscan_input\n\n'
            f'{exec_cmds} -f {output_path}/shell/lncRNA_prediction.sh -n 4\n\n'
            f'{plot_venn}\n\n'
            f'{lncRNA_results}\n\n'
            f'{create_lncRNA_gtf}\n\n'
            f'{lncRNA_exon_count}\n\n'
            f'{create_target_gtf}\n\n'
            f'{lncRNA_exp}\n\n'
            f'{target_exp}\n\n'
            f'{python} {output_path}/shell/lncRNA_target_prediction.py\n\n'
            f'{joint}\n\n'
            f'bash {output_path}/shell/lncRNA_classification.sh\n\n'
            f'{run_DE_enrich}'
        )
        o.write(cmd)
    system(f'chmod 755 {output_path}/shell/All_step.sh')
    click.echo(f'\033[32mCommands created successfully, please run "bash {output_path}/shell/All_step.sh".\033[0m', err=True)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('config', nargs=1, metavar='<yaml file|stdin>', type=click.File('r'), required=True)
@click.option('-V', '--version', 'version', help='Show author and version information.',
              is_flag=True, is_eager=True, expose_value=False, callback=displayer.version_info)
def run(config):
    """Strand specific (dUTP constructed library) pair end RNA-seq analysis pipeline (including lncRNA and its target gene prediction)."""
    main(config=config)


if __name__ == '__main__':
    run()
