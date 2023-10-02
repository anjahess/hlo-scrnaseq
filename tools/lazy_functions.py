"""

Collection of functions for any kind of repository.

@author: Anja Hess
@date: 2023-MAY-01

"""

import pandas as pd

def remove_duplicates(input_list):
    cleared_list = []
    for elem in input_list:
        if elem not in cleared_list:
            cleared_list.append(elem)
    return cleared_list


def adjust_list(input_list, list_valid):
    new_list = [e for e in input_list
                if e in list_valid]
    return new_list


def make_ensg2symbol(data_path):
    ensg2symbol = {}

    sources = ["features.tsv", 'mart_export.txt']

    for elem in sources:
        if elem.split(".")[1] == "tsv":
            input_df = pd.read_csv(
                f"{data_path}/{elem}",
                index_col=False,
                sep='\t')
        else:
            input_df = pd.read_csv(
                f"{data_path}/{elem}",
                index_col=False)

        input_df = input_df.to_dict(orient='index')

        for elem in input_df:
            if input_df[elem]['Gene name'] \
                    not in ensg2symbol:
                ensg2symbol.update({
                    input_df[elem]['Gene name']:
                        input_df[elem]['Gene stable ID']})
            try:
                if input_df[elem]['Gene name'] != \
                        input_df[elem]['Gene Synonym']:
                    ensg2symbol.update({
                        input_df[elem]['Gene Synonym']:
                            input_df[elem]['Gene stable ID']})
            except:
                continue
    return ensg2symbol


def make_feature_df(list_of_ensg, data_path):
    ens2symbol = make_ensg2symbol(data_path)

    for elem in list_of_ensg:
        if elem.split('.')[0] not in list(ens2symbol.keys()):
            ens2symbol.update({elem.split('.')[0]: "NO_SYMBOL"})

    feature_df = pd.DataFrame({'symbols':
                                   [str(ens2symbol[
                                            ensg.split('.')[0]])
                                    for ensg in list_of_ensg],
                               'gene_ids': list_of_ensg,
                               'feature_types':
                                   ["Gene Expression"
                                    for elem in list_of_ensg]
                               })
    feature_df = feature_df.set_index('symbols')
    return feature_df


def str2byte(list_of_strings):
    byte_list = []
    for elem in list_of_strings:
        elem = bytes(elem, 'utf-8')
        byte_list.append(elem)
    return byte_list


def byte2str(list_of_bytes):
    str_list = []
    for elem in list_of_bytes:
        print(elem)
        elem = elem.decode("utf-8")
        str_list.append(elem)
    return str_list


def logbook(runtitle="my_warning"):
    """
    :param runtitle: str
    :return: creates warnings > txt
    """
    import logging

    logging.basicConfig(
        filename=f"{runtitle}.txt",
        level=logging.INFO)
    logging.captureWarnings(True)

# END OF SCRIPT
