import pandas as pd
import py4cytoscape as p4c
from time import sleep
from py4cytoscape import CyError
pd.options.mode.chained_assignment = None
p4c.cytoscape_version_info()


def get_deseq_results(exp_filepath, genes_to_keep):
    experiment_de_df = pd.read_csv(exp_filepath, sep="\t")
    genes_to_keep_df = experiment_de_df[experiment_de_df['gene_id'].isin(genes_to_keep)]
    filtered_df_tmp_1 = experiment_de_df[experiment_de_df['log2FoldChange'].abs().gt(1)]
    filtered_df_tmp_2 = filtered_df_tmp_1[filtered_df_tmp_1['padj'].lt(0.05)]
    filtered_df = pd.concat([filtered_df_tmp_2, genes_to_keep_df]).drop_duplicates().reset_index(drop=True)
    return filtered_df

def get_network_name():
    network_suid = p4c.get_network_suid()
    network_name = p4c.get_network_name(network_suid)
    return network_name


def create_network(center_gene, specie, exp_name, exp_filepath, genes_to_keep):
    available_networks = p4c.get_network_list()
    main_network_name = center_gene + "_" + exp_name
    # if main_network_name in available_networks:
    #     p4c.set_current_network(main_network_name)
    #     sleep(0.5)
    # else:
    string_cmd_list = ['string protein query', 'query="', center_gene, '"', 'species=%s' % specie, 'cutoff=0.4',
                       "newNetName=%s" % main_network_name, "limit=1000"]
    string_cmd = " ".join(string_cmd_list)
    p4c.commands.commands_run(string_cmd)
    # p4c.set_current_network(main_network_name)
    sleep(0.5)
    string_cmd = 'string retrieve enrichment allNetSpecies="%s"'%specie
    p4c.commands.commands_run(string_cmd)
    sleep(0.5)
        # p4c.set_visual_style("de_genes_style")
    # mapped_cols = p4c.map_table_column('display name', specie, 'HGNC', 'Entrez Gene')
    # subnetwork_name = center_gene + "_" + exp_name
    # p4c.select_nodes(by_col= "display name" ,nodes=center_gene)
    # sleep(0.5)
    # p4c.select_first_neighbors()
    # p4c.create_subnetwork(nodes="selected", subnetwork_name=subnetwork_name)
    # sleep(0.5)
    deseq_df = get_deseq_results(exp_filepath, genes_to_keep)
    p4c.load_table_data(deseq_df, data_key_column='gene_id', table_key_column='display name')
    # p4c.export_
    string_cmd = 'string retrieve enrichment allNetSpecies="%s"'%specie
    p4c.commands.commands_run(string_cmd)
    sleep(0.5)
    # string_net = p4c.get_network_suid()
    # p4c.create_column_filter('log2filter', 'log2FoldChange', "", "IS_NOT", network=string_net)
    # p4c.apply_filter('log2filter', hide=True, network=string_net)
    p4c.set_visual_style("de_genes_style")
    return


def get_tables():
    exp_to_table = {}
    available_networks = p4c.get_network_list()
    for network in available_networks:
        # if "STRING" in network:
        #     continue
        exp = network.split("_")[1] + "_" + network.split("_")[2]
        gene = network.split("-")[1].strip().split("_")[0]
        if exp not in exp_to_table:
            exp_to_table[exp] = {}
        p4c.set_current_network(network)
        sleep(0.2)
        df = p4c.get_table_columns(columns=["gene_id", "log2FoldChange", ])
        df.dropna(subset=["log2FoldChange"], inplace=True)
        exp_to_table[exp][gene] = df["gene_id"].tolist()
    return exp_to_table

def exp_intersection(gene_bool_dict):
    exp_to_table = get_tables()
    intersection_list = []
    for exp in exp_to_table:
        for gene in gene_bool_dict:
            if gene_bool_dict[gene]:
                inclusion_list = exp_to_table[exp][gene]
                intersection_list.append(inclusion_list)
        genes_intersection = set.intersection(*map(set, intersection_list))
        for gene in gene_bool_dict:
            if not gene_bool_dict[gene]:
                exclusion_list = exp_to_table[exp][gene]
                final_genes_intersection = list(set(genes_intersection) - set(exclusion_list))
    return

def loop_experiments(experiments, specie, exp_filepath_variant, genes_to_keep):
    for center_gene in genes_to_keep:
        for exp_name in experiments:
            exp_filepath = exp_filepath_variant%(exp_name, exp_name)
            create_network(center_gene, specie, exp_name, exp_filepath, genes_to_keep)

if __name__ == '__main__':
    # experiments = ["r_nr", "nr_nrt", "r_rt", "rt_nrt"]
    # genes_to_keep = ['Braf', 'Akt1', 'Kras', 'Src']
    # exp_filepath_variant = "/home/direnc/results/tim/rnaseq_mice/%s/differential_abundance/tables/differential/condition_%s.deseq2.results.tsv"
    # specie = "Mus musculus"
    # loop_experiments(experiments, specie, exp_filepath_variant, genes_to_keep)
    tables_folder = "/home/direnc/results/tim/rnaseq_mice/gene_centered_results"
    gene_bool_dict = {"Kras" : True , "Akt1" : True , "Src" : True , "Braf" : False}
    exp_intersection(gene_bool_dict)