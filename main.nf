params.output = "./output/"
params.algorithm = "base" //"extra" or "all"
params.min_size = 20
params.timepoint = ""
params.meta_file = "./example/metadata.json"
params.count_files = "./example/data/1/quantification_table.tsv,./example/data/2/quantification_table.tsv,./example/data/3/quantification_table.tsv,./example/data/4/quantification_table.tsv,./example/data/5/quantification_table.tsv,./example/data/6/quantification_table.tsv,./example/data/7/quantification_table.tsv,./example/data/8/quantification_table.tsv,./example/data/10/quantification_table.tsv"

metadata = Channel.fromPath(params.meta_file)
file_list = params.count_files.tokenize(",")
file_channels = Channel.fromPath(file_list).collect()

// scripts
mosbi_script = Channel.fromPath("${projectDir}/universal_mosbi/universal.r")
protein_mapping = Channel.fromPath("${projectDir}/geneprot_mapping.csv")
join_table = Channel.fromPath("${projectDir}/semares_preprocessing/join_table.py")
metadata2table = Channel.fromPath("${projectDir}/semares_preprocessing/metadata2table.py")

// config files
data_config = Channel.fromPath("${projectDir}/config/data_table_config.json")
meta_data_config = Channel.fromPath("${projectDir}/config/meta_table_config.json")

process file_join {
    container "dockergenevention/pandas" // use docker conatainer
    publishDir params.output, mode: "copy"

    input:
    path script
    path metadata
    path config
    path file_channels, stageAs: "*/*"
    val file_path

    output:
    path "output/data.tsv"

    """
    python $script -m $metadata -c $config -f $file_channels -p $file_path
    """
}

process metadata_join {
    container "dockergenevention/pandas" // use docker conatainer
    publishDir params.output, mode: "copy"

    input:
    path script
    path metadata
    path config
    val file_path

    output:
    path "output/metadata.tsv"

    """
    python $script -m $metadata -c $config -p $file_path
    """
}

process mosbi {
    container "dockergenevention/mosbi" // use docker conatainer
    memory "8 GB"
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path count_file
    path meta_file
    path protein_mapping

    output:
    path "community*"
    path "Rplots.pdf"

    """
    Rscript $script_file ./ $count_file $meta_file ${params.algorithm} ${params.min_size} $protein_mapping ${params.timepoint}
    """
}

workflow {
  file_join(join_table, metadata, data_config, file_channels, params.count_files)
  metadata_join(metadata2table, metadata, meta_data_config, params.count_files)
  mosbi(mosbi_script, file_join.out, metadata_join.out, protein_mapping)
}
