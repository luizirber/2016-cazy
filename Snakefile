from functools import partial

include: "conf.py"

rule all:
    input:
        expand("outputs/common/{db}.x.{sample}-{depth}D.csv",
               db=DBS,
               sample=[1, 2, 3],
               depth=[2, 6, 15])

rule download_dbcan:
    output:
        "inputs/db/cazy.dbcan.fa"
    shell: """
        wget http://csbl.bmb.uga.edu/dbCAN/download/CAZyDB.07152016.fa -O {output}
    """

rule prepare_db_diamond_seqs:
    input: "outputs/merged/{sample}-{depth}D.assembled.fasta.gz"
    output: "outputs/diamond/db/samples/{sample}-{depth}D.dmnd"
    shell: """
        diamond makedb --in {input} -d {output}
    """

rule prepare_db_diamond_cazy:
    input: "inputs/db/cazy.{dbid}.fa"
    output: "outputs/diamond/db/cazy/{dbid}.dmnd"
    shell: """
        diamond makedb --in {input} -d {output}
    """

rule diamond_search_seqs:
    input:
       db="outputs/diamond/db/samples/{sample}-{depth}D.dmnd",
       query="inputs/db/cazy.{dbid}.fa"
    output: "outputs/diamond/{sample}-{depth}D/{dbid}.daa"
    threads: 4
    shell: """
        diamond blastp -k 100 --sensitive \
                       -e {EVALUE} \
                       -d {input.db} -q {input.query} \
                       -a {output} -p {threads}
    """

rule diamond_search_db:
    input:
       db="outputs/diamond/db/cazy/{dbid}.dmnd",
       query="outputs/merged/{sample}-{depth}D.assembled.fasta.gz"
    output: "outputs/diamond/{dbid}/{sample}-{depth}D.daa"
    threads: 4
    shell: """
        diamond blastx -k 100 --sensitive \
                       -e {EVALUE} \
                       -d {input.db} -q {input.query} \
                       -a {output} -p {threads}
    """

rule diamond_convert_output:
    input:
      "outputs/diamond/{db}/{query}.daa"
    output:
      "outputs/diamond/{db}/{query}.out.gz"
    shell:
        "diamond view -a {input} | gzip > {output}"

rule find_common_occurences:
    input:
      "outputs/diamond/{db}/{query}.out.gz"
    output:
      "outputs/common/{db}.x.{query}.csv"
    shell:
      "touch {output}"

rule download_pear:
    output: "sw/pear"
    shell: """
        mkdir -p sw
        curl -O http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz
        tar xf pear-0.9.10-bin-64.tar.gz
        mv pear-0.9.10-bin-64/pear-0.9.10-bin-64 {output}
        rm -rf pear-0.9.10-bin-64*
    """

def raw_inputs(w, direction=1):
    return "inputs/Metagenomes_Harbour_Fagans_BAY/{sample}-{depth}D_C3ETWACXX_{barcode}_L007_R{d}.fastq.gz".format(
          barcode=BARCODE[w.sample + '-' + w.depth], d=direction, **w)

rule merge_pairs:
    input:
        forward=partial(raw_inputs, direction=1),
        reverse=partial(raw_inputs, direction=2)
    output: "outputs/merged/{sample}-{depth}D.assembled.fastq"
    params: outprefix="outputs/merged/{sample}-{depth}D"
    threads: 4
    shell: """
        sw/pear -f {input.forward} -r {input.reverse} -o {PRJ_ROOT}/{params.outprefix} \
                 -p 0.05 -q 25 -y 5G -j {threads} -v 4 -g 2 -e
    """

rule fastq_to_fasta:
    input: "outputs/merged/{sample}-{depth}D.assembled.fastq"
    output: "outputs/merged/{sample}-{depth}D.assembled.fasta.gz"
    shell: "fastq-to-fasta.py -n --gzip -o {output} {input}"

ruleorder: prepare_db_diamond_seqs > prepare_db_diamond_cazy
