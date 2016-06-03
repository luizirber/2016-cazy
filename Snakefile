from functools import partial


PRJ_ROOT = next(shell("readlink -e .", iterable=True))
SUPERFOCUS_DIR = "{0}/sw/SUPERFOCUS_0.24".format(PRJ_ROOT)


rule all:
    input:
        expand("outputs/diamond/{db}/{sample}-{depth}D.daa",
               db=["dbcan", "tools.aa", "tools.ce", #"tools.gh",
                   "tools.gt", "tools.pl"],
               sample=[1, 2, 3],
               depth=[2, 6, 15]
               ),
        expand("outputs/diamond/{sample}-{depth}D/{db}.daa",
               db=["dbcan", "tools.aa", "tools.ce", #"tools.gh",
                   "tools.gt", "tools.pl"],
               sample=[1, 2, 3],
               depth=[2, 6, 15]
               )

rule download_dbcan:
    output:
        "inputs/db/cazy.dbcan.fa"
    shell: """
        wget http://csbl.bmb.uga.edu/dbCAN/download/CAZyDB.03172015.fa -O {output}
    """

rule download_cazy_tools:
    output:
        "inputs/db/cazy.tools.{family}.fa"
    params: family="{family}"
    shell: """
        wget http://bit.ly/cazy-{params.family} -O {output}.zip
        unzip -p {output}.zip > {output}
        rm {output}.zip
    """

rule prepare_db_seqs:
    input: "outputs/merged/{sample}-{depth}D.assembled.fasta.gz"
    output: "outputs/diamond/db/samples/{sample}-{depth}D.dmnd"
    shell: """
        diamond makedb --in {input} -d {output}
    """

rule prepare_db_cazy:
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
        diamond blastp -d {input.db} -q {input.query} -a {output} -p {threads}
    """

rule diamond_search_db:
    input:
       db="outputs/diamond/db/cazy/{dbid}.dmnd",
       query="outputs/merged/{sample}-{depth}D.assembled.fasta.gz"
    output: "outputs/diamond/{dbid}/{sample}-{depth}D.daa"
    threads: 4
    shell: """
        diamond blastx -d {input.db} -q {input.query} -a {output} -p {threads}
    """


rule download_pear:
    output: "sw/pear"
    shell: """
        mkdir -p sw
        git clone https://github.com/xflouris/PEAR.git
        cp PEAR/bin/pear-0.9.6-bin-64 {output}
        rm -rf PEAR
    """

BARCODE = {
    "1-2":  "CAGATC",
    "1-6":  "AGTCAA",
    "1-15": "TAGCTT",
    "2-2":  "ACTTGA",
    "2-6":  "AGTTCC",
    "2-15": "GGCTAC",
    "3-2":  "GATCAG",
    "3-6":  "ATGTCA",
    "3-15": "CTTGTA",
}


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

rule  fastq_to_fasta:
    input: "outputs/merged/{sample}-{depth}D.assembled.fastq"
    output: "outputs/merged/{sample}-{depth}D.assembled.fasta.gz"
    shell: "fastq-to-fasta.py -n --gzip -o {output} {input}"

rule download_SUPERFOCUS:
    output: "sw/SUPERFOCUS_0.24/superfocus.py"
    shell: """
        mkdir sw
        wget -c http://downloads.sourceforge.net/project/superfocus/standalone/SUPERFOCUS_0.24.zip -P sw
        cd sw && unzip SUPERFOCUS_0.24.zip
        cd {SUPERFOCUS_DIR} && python superfocus__downloadDB.py rapsearch blast diamond
    """.format(SUPERFOCUS_DIR=SUPERFOCUS_DIR)

ruleorder: prepare_db_seqs > prepare_db_cazy
