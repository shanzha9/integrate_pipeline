import csv
from collections import Counter
import os


def locate_reply_elements_in_list(input_list, single_result=True):
    old_list_copy = []
    reply_index = {}

    for index, i in enumerate(input_list):
        if i in old_list_copy:
            index_list = reply_index.get(i)
            if index_list is None:
                pass
            else:
                index_list.append(index)
                reply_index.update({i: index_list})
        else:
            reply_index.update({i: [index]})
            old_list_copy.append(i)

    if single_result:
        reply_index_result = {}
        for k, v in reply_index.items():
            if len(v) == 1:
                pass
            else:
                reply_index_result.setdefault(k, v)
        return reply_index_result
    else:
        return reply_index

def rwFiles(fileName='result.csv', model='r', line=None):
    """
    Read and write files like CSV
    :param fileName: input file name
    :param model: r or a+
    :param line: write lines
    :return: List data
    """
    if model == 'r':
        with open(fileName, mode=model) as _f:
            _data = _f.read().split('\n')
        if len(_data[-1]) == 0:
            _data.pop(-1)
        return _data
    else:
        with open(fileName, mode=model) as _f:
            line = line + '\n'
            _f.write(line)


def count_umi(input_file_path=None, sample_name=None, output_file_name=None):
    
    line_list = []
    barcode_list = []
    with open(input_file_path) as f:
        raw_data = csv.reader(f)
        for index, line in enumerate(raw_data):
            if index == 0:
                index_barcode, index_hc, index_chain, index_prod, index_reads, index_umis, index_vgene, index_dgene, index_jgene, index_cgene, index_cdr3, index_cdr3nt, index_raw_clonotypeid = line.index(
                    "barcode"), line.index("high_confidence"), line.index("chain"), line.index(
                    "productive"), line.index("reads"), line.index("umis"), line.index("v_gene"), line.index("d_gene"), line.index("j_gene"), line.index("c_gene"), line.index("cdr3"), line.index("cdr3_nt"), line.index("raw_clonotype_id")
                if "BCR" in input_file_path:
                    line = ['sample_name', 'bcr_umi', 'bcr_vgene', 'bcr_dgene', 'bcr_jgene', 'bcr_cgene', 'bcr_cdr3', 'bcr_cdr3nt', 'bcr_rawClonotypeid']
                else:
                    line = ['sample_name', 'tcr_umi', 'tcr_vgene', 'tcr_dgene', 'tcr_jgene', 'tcr_cgene', 'tcr_cdr3', 'tcr_cdr3nt', 'tcr_rawClonotypeid']
                line = ','.join(line)
                rwFiles(output_file_name, model='a+', line=line)
            else:
                if line[index_hc].upper() == "TRUE":
                    if line[index_prod].upper() == "TRUE":
                        if line[index_chain].upper() != "MULTI":
                            line_list.append(line)
                            barcode_list.append(line[index_barcode])
    
    barcode_copy_index = locate_reply_elements_in_list(barcode_list)
    length_ = []
    for _, v in barcode_copy_index.items():
        length_.append(len(v))
    
    deal_list = []
    for i in line_list:
        if i[index_barcode] in barcode_copy_index.keys():
            """重复barcode处理:
            1.如果重复的barcode数量为2，直接将UMI相加
            2.如果重复的barcode数量大于2，需确定重复的是哪一类？再取reads最大对应的类，与另一类相加
            """
            if len(barcode_copy_index.get(i[index_barcode])) == 2:
                if i[index_barcode] in deal_list:
                    pass
                else:
                    [index_x, index_y] = barcode_copy_index.get(i[index_barcode])
                    umi = int(line_list[index_x][index_umis]) + int(line_list[index_y][index_umis])
                    deal_list.append(i[index_barcode])
                    _name = '{}_{}'.format(sample_name, i[index_barcode])
                    umi = str(umi)
                    v_gene = i[index_vgene]
                    d_gene = i[index_dgene]
                    j_gene = i[index_jgene]
                    c_gene = i[index_cgene]
                    cdr3 = i[index_cdr3]
                    cdr3nt = i[index_cdr3nt]
                    raw_clonotypeid = i[index_raw_clonotypeid]
                    line = [_name, umi, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3nt, raw_clonotypeid]
                    line = ','.join(line)
                    rwFiles(output_file_name, model='a+', line=line)
            else:
                if i[index_barcode] in deal_list:
                    pass
                else:
                    # 把重复的chain放入新的list寻找重复的元素
                    reply_site_list = []
                    wait_deal_line_list = []
                    for j in barcode_copy_index.get(i[index_barcode]):
                        reply_site_list.append(line_list[j][index_chain])
                        wait_deal_line_list.append(line_list[j])
                    chain_copy_index = locate_reply_elements_in_list(input_list=reply_site_list,
                                                                           single_result=False)
                    umi = 0
                    for _, v in chain_copy_index.items():
                        max_reads = 0
                        max_umi = 0
                        for k in v:
                            reads = int(wait_deal_line_list[k][index_reads])
                            umi_ = int(wait_deal_line_list[k][index_umis])
                            if reads > max_reads:
                                max_reads = reads
                                max_umi = umi_
                        umi = umi + max_umi
                    _name = '{}_{}'.format(sample_name, i[index_barcode])
                    umi = str(umi)
                    v_gene = i[index_vgene]
                    d_gene = i[index_dgene]
                    j_gene = i[index_jgene]
                    c_gene = i[index_cgene]
                    cdr3 = i[index_cdr3]
                    cdr3nt = i[index_cdr3nt]
                    raw_clonotypeid = i[index_raw_clonotypeid]
                    line = [_name, umi, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3nt, raw_clonotypeid]
                    line = ','.join(line)
                    rwFiles(output_file_name, model='a+', line=line)
                    deal_list.append(i[index_barcode])
        else:
            # 不做处理
            _name = '{}_{}'.format(sample_name, i[index_barcode])
            umi = str(i[index_umis])
            v_gene = i[index_vgene]
            d_gene = i[index_dgene]
            j_gene = i[index_jgene]
            c_gene = i[index_cgene]
            cdr3 = i[index_cdr3]
            cdr3nt = i[index_cdr3nt]
            raw_clonotypeid = i[index_raw_clonotypeid]
            line = [_name, umi, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3nt, raw_clonotypeid]
            line = ','.join(line)
            rwFiles(output_file_name, model='a+', line=line)
