from datetime import date


##### Configuration #####
config["reference_genbank"] = config["reference_genbank"].rstrip("/")
config["fasta"] = config["fasta"].rstrip("/")
config["output_path"] = config["output_path"].rstrip("/")

today = date.today().strftime("%m%d%y")

##### Workflow #####
rule all:
    input:
        expand(config["output_path"] + "/{date}.fasta",date=today)

rule parse_reference:
    input:
        config["reference_genbank"]
    output:
        config["reference_genbank"] + ".json"
    log:
        config["output_path"] + "/refparser.log"
    shell:
        genofunk refparser -g {input} -o {output} -v &> {log}

rule annotate:
    input:
        reference_genbank = rules.parse_reference.output,
        consensus_fasta = config["fasta"]
    output:
        config["fasta"] + ".coordinates",
        config["fasta"] + ".edits",
    log:
        config["output_path"] + "/annotate.log"
    shell:
        genofunk annotate -g {input} -o {output} -v &> {log}