#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MEME XML解析与基因二进制编码生成脚本

功能：
1. 解析MEME.xml文件
2. 提取每个基因的motif信息  
3. 为每个基因生成二进制编码（从左往右依次表示第几个Motif的存在与否）
4. 1表示存在，0表示不存在
"""

import xml.etree.ElementTree as ET
import argparse
import sys
import os
from collections import defaultdict


def parse_meme_xml(xml_file):
    """
    解析MEME XML文件
    
    Args:
        xml_file: MEME输出的XML文件路径
        
    Returns:
        tuple: (motifs_info, gene_motifs, sequence_info)
        - motifs_info: 字典，键为motif_id，值为motif详细信息
        - gene_motifs: 字典，键为基因名，值为该基因包含的motif集合
        - sequence_info: 字典，键为sequence_id，值为序列信息
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        motifs_info = {}
        gene_motifs = defaultdict(set)
        sequence_info = {}
        
        print(f"正在解析文件: {xml_file}")
        
        # 首先解析序列信息
        training_set = root.find('training_set')
        if training_set is not None:
            sequences = training_set.findall('sequence')
            print(f"发现 {len(sequences)} 个序列")
            
            for seq in sequences:
                seq_id = seq.get('id')
                seq_name = seq.get('name')
                seq_length = seq.get('length')
                seq_weight = seq.get('weight')
                
                sequence_info[seq_id] = {
                    'id': seq_id,
                    'name': seq_name,
                    'length': seq_length,
                    'weight': seq_weight
                }
                
                # 初始化每个基因的motif集合
                gene_motifs[seq_name] = set()
        else:
            print("警告: 未找到training_set元素")
        
        # 解析motif信息
        motifs_element = root.find('motifs')
        if motifs_element is not None:
            motifs = motifs_element.findall('motif')
            print(f"发现 {len(motifs)} 个motif")
            
            for motif in motifs:
                motif_id = motif.get('id')
                motif_name = motif.get('name') 
                width = motif.get('width')
                sites = motif.get('sites')
                ic = motif.get('ic')
                re = motif.get('re')
                llr = motif.get('llr')
                p_value = motif.get('p_value')
                e_value = motif.get('e_value')
                
                motifs_info[motif_id] = {
                    'id': motif_id,
                    'name': motif_name,
                    'width': width,
                    'sites': sites,
                    'ic': ic,
                    're': re,
                    'llr': llr,
                    'p_value': p_value,
                    'e_value': e_value
                }
                
                # 解析contributing_sites来确定哪些基因包含此motif
                contributing_sites_elem = motif.find('contributing_sites')
                if contributing_sites_elem is not None:
                    contributing_sites = contributing_sites_elem.findall('contributing_site')
                    
                    for site in contributing_sites:
                        sequence_id = site.get('sequence_id')
                        position = site.get('position')
                        strand = site.get('strand')
                        pvalue = site.get('pvalue')
                        
                        if sequence_id in sequence_info:
                            gene_name = sequence_info[sequence_id]['name']
                            gene_motifs[gene_name].add(motif_id)
                        else:
                            print(f"警告: 序列ID '{sequence_id}' 在training_set中未找到")
        else:
            print("错误: 未找到motifs元素")
            sys.exit(1)
        
        print(f"成功解析 {len(sequence_info)} 个序列和 {len(motifs_info)} 个motif")
        return motifs_info, gene_motifs, sequence_info
        
    except ET.ParseError as e:
        print(f"XML解析错误: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"文件读取错误: {e}")
        sys.exit(1)


def generate_binary_encoding(motifs_info, gene_motifs):
    """
    为每个基因生成二进制编码
    
    Args:
        motifs_info: motif信息字典
        gene_motifs: 基因-motif对应关系字典
        
    Returns:
        tuple: (gene_encodings, sorted_motif_ids)
        - gene_encodings: 基因名到二进制编码的映射
        - sorted_motif_ids: 排序的motif ID列表
    """
    # 按motif ID中的数字排序，确保motif_1, motif_2, motif_10的顺序正确
    def motif_sort_key(motif_id):
        """提取motif ID中的数字进行排序"""
        try:
            if motif_id.startswith('motif_'):
                return int(motif_id.split('_')[1])
            else:
                return float('inf')  # 将非标准格式放到最后
        except (IndexError, ValueError):
            return float('inf')
    
    sorted_motif_ids = sorted(motifs_info.keys(), key=motif_sort_key)
    total_motifs = len(sorted_motif_ids)
    
    print(f"\nMotif排序列表 (共{total_motifs}个):")
    for i, motif_id in enumerate(sorted_motif_ids, 1):
        motif_name = motifs_info[motif_id].get('name', 'N/A')
        sites_count = motifs_info[motif_id].get('sites', 'N/A')
        print(f"  位置 {i}: {motif_id} ({motif_name}) - {sites_count} sites")
    
    gene_encodings = {}
    
    # 确保所有基因都有编码，即使没有任何motif
    all_genes = set(gene_motifs.keys())
    
    for gene_name in all_genes:
        # 初始化二进制编码（全为0）
        binary_code = ['0'] * total_motifs
        
        # 检查每个motif是否在该基因中存在
        for i, motif_id in enumerate(sorted_motif_ids):
            if motif_id in gene_motifs[gene_name]:
                binary_code[i] = '1'
        
        # 转为字符串
        gene_encodings[gene_name] = ''.join(binary_code)
    
    return gene_encodings, sorted_motif_ids


def write_results(gene_encodings, sorted_motif_ids, motifs_info, output_file):
    """
    将结果写入输出文件
    
    Args:
        gene_encodings: 基因编码字典
        sorted_motif_ids: 排序的motif ID列表
        motifs_info: motif详细信息
        output_file: 输出文件路径
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入表头
        f.write("# MEME基因二进制编码结果\n")
        f.write(f"# 分析时间: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# 总共 {len(sorted_motif_ids)} 个Motif, {len(gene_encodings)} 个基因\n")
        f.write("# 编码格式：从左到右依次表示第1到第N个Motif的存在与否\n")
        f.write("# 1=存在，0=不存在\n\n")
        
        # 写入motif对应位置说明
        f.write("# Motif位置对应表:\n")
        for i, motif_id in enumerate(sorted_motif_ids, 1):
            motif_name = motifs_info[motif_id].get('name', 'N/A')
            width = motifs_info[motif_id].get('width', 'N/A')
            sites = motifs_info[motif_id].get('sites', 'N/A')
            e_value = motifs_info[motif_id].get('e_value', 'N/A')
            f.write(f"# 位置 {i}: {motif_id} ({motif_name}) - 宽度:{width}bp, 位点数:{sites}, E-value:{e_value}\n")
        f.write("\n")
        
        # 写入列标题
        f.write("基因名\t二进制编码\t包含的Motif数\t包含的Motif列表\n")
        
        # 写入每个基因的编码结果
        for gene_name in sorted(gene_encodings.keys()):
            encoding = gene_encodings[gene_name]
            
            # 找出该基因包含的motif
            contained_motifs = []
            motif_count = 0
            for i, bit in enumerate(encoding):
                if bit == '1':
                    contained_motifs.append(sorted_motif_ids[i])
                    motif_count += 1
            
            motif_list = ','.join(contained_motifs) if contained_motifs else 'None'
            f.write(f"{gene_name}\t{encoding}\t{motif_count}\t{motif_list}\n")


def print_summary(gene_encodings, sorted_motif_ids, motifs_info):
    """
    打印结果摘要
    
    Args:
        gene_encodings: 基因编码字典
        sorted_motif_ids: 排序的motif ID列表
        motifs_info: motif详细信息
    """
    print("\n" + "="*80)
    print("结果摘要")
    print("="*80)
    print(f"总Motif数量: {len(sorted_motif_ids)}")
    print(f"总基因数量: {len(gene_encodings)}")
    
    # 统计每个motif在多少个基因中出现
    motif_counts = defaultdict(int)
    for gene_name, encoding in gene_encodings.items():
        for i, bit in enumerate(encoding):
            if bit == '1':
                motif_counts[sorted_motif_ids[i]] += 1
    
    print(f"\n各Motif出现频率:")
    for motif_id in sorted_motif_ids:
        count = motif_counts[motif_id]
        percentage = (count / len(gene_encodings)) * 100 if len(gene_encodings) > 0 else 0
        motif_name = motifs_info[motif_id].get('name', 'N/A')
        e_value = motifs_info[motif_id].get('e_value', 'N/A')
        print(f"  {motif_id} ({motif_name}): {count}/{len(gene_encodings)}个基因 ({percentage:.1f}%) - E-value: {e_value}")
    
    # 统计基因中包含motif的分布
    motif_per_gene = [sum(1 for bit in encoding if bit == '1') for encoding in gene_encodings.values()]
    if motif_per_gene:
        avg_motifs = sum(motif_per_gene) / len(motif_per_gene)
        max_motifs = max(motif_per_gene)
        min_motifs = min(motif_per_gene)
        
        print(f"\n基因中Motif分布:")
        print(f"  平均每个基因包含Motif数: {avg_motifs:.2f}")
        print(f"  最多包含Motif数: {max_motifs}")
        print(f"  最少包含Motif数: {min_motifs}")
    
    # 显示前10个基因的编码示例
    print(f"\n前10个基因编码示例:")
    for i, (gene_name, encoding) in enumerate(sorted(gene_encodings.items())[:10]):
        motif_count = sum(1 for bit in encoding if bit == '1')
        print(f"  {gene_name}: {encoding} (包含{motif_count}个motif)")
    
    if len(gene_encodings) > 10:
        print(f"  ... (还有{len(gene_encodings)-10}个基因)")


def validate_xml_structure(xml_file):
    """
    验证XML文件结构的完整性
    
    Args:
        xml_file: XML文件路径
        
    Returns:
        bool: 是否为有效的MEME XML文件
    """
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        if root.tag != 'MEME':
            print("错误: 这不是一个有效的MEME XML文件")
            return False
        
        training_set = root.find('training_set')
        motifs = root.find('motifs')
        
        if training_set is None:
            print("错误: XML文件中缺少training_set元素")
            return False
            
        if motifs is None:
            print("错误: XML文件中缺少motifs元素")
            return False
        
        sequences = training_set.findall('sequence')
        motif_elements = motifs.findall('motif')
        
        if len(sequences) == 0:
            print("错误: 未找到任何序列信息")
            return False
            
        if len(motif_elements) == 0:
            print("错误: 未找到任何motif信息")
            return False
        
        print(f"XML结构验证通过: 发现{len(sequences)}个序列, {len(motif_elements)}个motif")
        return True
        
    except Exception as e:
        print(f"XML结构验证失败: {e}")
        return False


def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(
        description="解析MEME XML文件并生成基因的motif二进制编码",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python meme_binary_encoder.py -i meme.xml -o motif_binary_encoding.txt
  python meme_binary_encoder.py -i meme.xml -o results.txt --verbose
  python meme_binary_encoder.py -i meme.xml -o results.txt --validate
        """
    )
    
    parser.add_argument('-i', '--input', 
                       required=True,
                       help='输入的MEME XML文件路径')
    
    parser.add_argument('-o', '--output', 
                       required=True,
                       help='输出文件路径')
    
    parser.add_argument('--verbose', 
                       action='store_true',
                       help='显示详细信息')
    
    parser.add_argument('--validate', 
                       action='store_true',
                       help='验证XML文件结构')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 '{args.input}' 不存在!")
        sys.exit(1)
    
    # 验证XML结构
    if args.validate:
        if not validate_xml_structure(args.input):
            sys.exit(1)
    
    # 解析XML文件
    motifs_info, gene_motifs, sequence_info = parse_meme_xml(args.input)
    
    if not motifs_info:
        print("错误: 未从XML文件中找到任何motif信息!")
        sys.exit(1)
    
    if not gene_motifs:
        print("错误: 未从XML文件中找到任何基因-motif关联信息!")
        sys.exit(1)
    
    # 生成二进制编码
    gene_encodings, sorted_motif_ids = generate_binary_encoding(motifs_info, gene_motifs)
    
    # 写入结果文件
    write_results(gene_encodings, sorted_motif_ids, motifs_info, args.output)
    
    # 打印摘要
    if args.verbose:
        print_summary(gene_encodings, sorted_motif_ids, motifs_info)
    
    print(f"\n结果已保存到: {args.output}")
    print("="*50)
    print("编码说明:")
    print("- 每行代表一个基因")
    print("- 二进制编码从左到右对应motif_1, motif_2, ...")
    print("- 1表示基因包含该motif，0表示不包含")
    print("- 详细的motif信息请查看输出文件头部")
    print("="*50)
    print("完成!")


if __name__ == "__main__":
    main()