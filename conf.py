PRJ_ROOT = next(shell("readlink -e .", iterable=True))
DBS = ["dbcan"]

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

EVALUE = "1e-10"
