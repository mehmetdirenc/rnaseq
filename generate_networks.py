import pandas as pd
import py4cytoscape as p4c
from time import sleep
import os
pd.options.mode.chained_assignment = None
p4c.cytoscape_version_info()
from itertools import combinations
from pyvenn import *


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
    if "STRING network - %s"%main_network_name in available_networks:
        return
    string_cmd_list = ['string protein query', 'query="', center_gene, '"', 'species=%s' % specie, 'cutoff=0.4',
                       "newNetName=%s" % main_network_name, "limit=1000"]
    string_cmd = " ".join(string_cmd_list)
    p4c.commands.commands_run(string_cmd)
    sleep(0.5)
    string_cmd = 'string retrieve enrichment allNetSpecies="%s"'%specie
    p4c.commands.commands_run(string_cmd)
    sleep(0.5)
    deseq_df = get_deseq_results(exp_filepath, genes_to_keep)
    p4c.load_table_data(deseq_df, data_key_column='gene_id', table_key_column='display name')
    string_cmd = 'string retrieve enrichment allNetSpecies="%s"'%specie
    p4c.commands.commands_run(string_cmd)
    sleep(0.5)
    # p4c.set_visual_style("de_genes_style")
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
        sleep(0.1)
        df = p4c.get_table_columns(columns=["gene_id", "log2FoldChange", ])
        df.dropna(subset=["log2FoldChange"], inplace=True)
        exp_to_table[exp][gene] = df["gene_id"].tolist()
    return exp_to_table

def exp_intersection(gene_bool_dict, main_res_folder):
    exp_to_table = get_tables()
    for exp in exp_to_table:
        intersection_list = []
        if exp == "rt_nrt":
            print("asd")
        create_venn(exp_to_table[exp], exp, main_res_folder)
        for gene in gene_bool_dict:
            if gene_bool_dict[gene]:
                inclusion_list = exp_to_table[exp][gene]
                intersection_list.append(inclusion_list)
        genes_intersection = set.intersection(*map(set, intersection_list))
        for gene in gene_bool_dict:
            if not gene_bool_dict[gene]:
                exclusion_list = exp_to_table[exp][gene]
                final_genes_intersection = list(set(genes_intersection) - set(exclusion_list))
        write_gene_list(final_genes_intersection, main_res_folder, exp)
    return

def exp_intersection_all_combinations(gene_bool_dict, main_res_folder, areas_of_interest):
    gene_set = [i for i in gene_bool_dict.keys() if gene_bool_dict[i] == True]
    exp_to_table = get_tables()

    for exp in exp_to_table:
        all_collections = create_venn(exp_to_table[exp], exp, main_res_folder)
        write_gene_list_all_collections(all_collections, main_res_folder, exp, areas_of_interest)
        # current_exp_dict = exp_to_table[exp]
        # current_exp_dict_set = {key: set(value) for key, value in current_exp_dict.items()}
        # for i in gene_bool_dict.keys():
        #     if gene_bool_dict[i] == False:
        #         current_exp_dict_set.pop(i)
        # intersections = {}
        # for i in range(1, len(current_exp_dict_set) + 1):
        #     for combo in combinations(current_exp_dict_set.keys(), i):
        #         intersected_genes = set.intersection(*(current_exp_dict_set[set_name] for set_name in combo))
        #         for gene in gene_bool_dict:
        #             if not gene_bool_dict[gene]:
        #                 exclusion_list = exp_to_table[exp][gene]
        #                 final_genes_intersection = list(intersected_genes - set(exclusion_list))
        #         intersections[combo] = final_genes_intersection
        #
        # write_gene_list(intersections, main_res_folder, exp)
    return




def create_venn(exp_dict, exp, main_res_folder):
    res_filepath = os.path.join(main_res_folder, "venn", exp)
    exp_dict_set = {}
    genes = []
    exp_dict = dict(sorted(exp_dict.items()))
    for gene in exp_dict:
        exp_dict_set[gene] = set(exp_dict[gene])
        genes.append(gene)
    labels1 = list(exp_dict_set.values())
    labels, all_collections = get_labels(labels1)
    fig, ax = venn4(labels, names=list(exp_dict_set.keys()))
    # plt.show()
    plt.savefig(res_filepath)
    plt.close()
    return all_collections

def loop_experiments(experiments, specie, exp_filepath_variant, genes_to_keep):
    for center_gene in genes_to_keep:
        for exp_name in experiments:
            exp_filepath = exp_filepath_variant%(exp_name, exp_name)
            create_network(center_gene, specie, exp_name, exp_filepath, genes_to_keep)

def write_gene_list(intersection_combo, main_res_folder, exp):
    result_filepath = os.path.join(main_res_folder, "gene_lists", exp + ".tsv")
    # with open(result_filepath, 'w') as rf:
    max_length = max(len(lst) for lst in intersection_combo.values())
    intersection_combo2 = {}
    # Pad lists with None to have the same length
    for key in intersection_combo:
        new_key = "_".join(key) + "_" + str(len(intersection_combo[key]))
        intersection_combo2[new_key] = intersection_combo[key] + [None] * (max_length - len(intersection_combo[key]))

    # Convert the dictionary to a DataFrame
    intersection_df = pd.DataFrame(intersection_combo2).fillna("")
    # intersection_df = pd.DataFrame(intersection_combo.items())
    # rf.write("#Gene\n")
    intersection_df.to_csv(result_filepath, index=False, sep="\t")
        # for i in intersection_combo:
        #     rf.write("%s\t"%i)
        # rf.write("\n")
        # for gene in gene_list:
        #     rf.write("%s\n" % gene)


def write_gene_list_all_collections(all_collections, main_res_folder, exp, areas_of_interest):
    result_filepath = os.path.join(main_res_folder, "gene_lists", exp + ".tsv")
    intersection = {}
    for area_name in areas_of_interest:
        all_genes = []
        for area in areas_of_interest[area_name]:
            for gene in all_collections[area]:
                all_genes.append(gene)
        intersection[area_name] = all_genes
    # with open(result_filepath, 'w') as rf:
    # max_length = max(len(lst) for lst in intersection_combo.values())
    # intersection_combo2 = {}
    # # Pad lists with None to have the same length
    # for key in intersection_combo:
    #     new_key = "_".join(key) + "_" + str(len(intersection_combo[key]))
    #     intersection_combo2[new_key] = intersection_combo[key] + [None] * (ma   x_length - len(intersection_combo[key]))

    # Convert the dictionary to a DataFrame
    pad_dict_list(intersection, "")
    intersection_df = pd.DataFrame(intersection)
    # intersection_df = pd.DataFrame(intersection_combo.items())
    # rf.write("#Gene\n")
    intersection_df.to_csv(result_filepath, index=False, sep="\t")
        # for i in intersection_combo:
        #     rf.write("%s\t"%i)
        # rf.write("\n")
        # for gene in gene_list:
        #     rf.write("%s\n" % gene)

def pad_dict_list(dict_list, padel):
    lmax = 0
    for lname in dict_list.keys():
        lmax = max(lmax, len(dict_list[lname]))
    for lname in dict_list.keys():
        ll = len(dict_list[lname])
        if  ll < lmax:
            dict_list[lname] += [padel] * (lmax - ll)
    return dict_list

if __name__ == '__main__':
    # experiments = ["nr_r", "nr_nrt", "r_rt", "nrt_rt"]
    # genes_to_keep = ['Braf', 'Akt1', 'Kras', 'Src']
    # exp_filepath_variant = "/home/direnc/results/tim/rnaseq_mice/%s/differential_abundance/tables/differential/condition_%s.deseq2.results.tsv"
    # specie = "Mus musculus"
    # loop_experiments(experiments, specie, exp_filepath_variant, genes_to_keep)
    main_res_folder = "/home/direnc/results/tim/rnaseq_mice/gene_centered_results"
    gene_bool_dict = {"Kras" : True , "Akt1" : True , "Src" : True , "Braf" : False}
    areas_of_interest = {"kras" : ["0010"], "kras_src" : ["0010", "0011"]}
    # exp_intersection(gene_bool_dict, main_res_folder)
    exp_intersection_all_combinations(gene_bool_dict, main_res_folder, areas_of_interest)