#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
拓扑编码和Motif编码合并脚本

功能：
1. 读取拓扑编码文件（topology_codes_summary.txt）
2. 读取motif二进制编码文件（gene_motif_binary.txt）
3. 按照拓扑文件中的顺序合并数据
4. 输出格式化结果文件
"""

import argparse
import sys
import os
from collections import defaultdict, OrderedDict


def parse_topology_file(file_path, verbose=False):
    """
    解析拓扑编码文件
    
    Args:
        file_path: 拓扑编码文件路径
        verbose: 是否显示详细信息
        
    Returns:
        dict: {基因名: OrderedDict{物种名: 拓扑编码}}
    """
    topology_data = defaultdict(OrderedDict)
    current_gene = None
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # 跳过空行和分隔线
                if not line or line.startswith('='):
                    continue
                
                # 检测基因名行
                if line.startswith('Gene:') or line.startswith('基因:'):
                    current_gene = line.split(':', 1)[1].strip()
                    # 处理包含特殊字符的基因名
                    if verbose:
                        print(f"   发现基因: {current_gene}")
                    continue
                
                # 跳过表头
                if '物种名称' in line and '拓扑编码' in line:
                    continue
                
                # 解析数据行
                if current_gene and '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        species_name = parts[0].strip()
                        topo_code = parts[1].strip()
                        topology_data[current_gene][species_name] = topo_code
                    
        print(f"成功解析拓扑文件: {file_path}")
        for gene, species_data in topology_data.items():
            print(f"   基因 {gene}: {len(species_data)} 个物种")
            
        return topology_data
        
    except FileNotFoundError:
        print(f"拓扑文件不存在: {file_path}")
        return {}
    except Exception as e:
        print(f"解析拓扑文件失败: {e}")
        return {}


def parse_motif_file(file_path):
    """
    解析motif二进制编码文件
    
    Args:
        file_path: motif编码文件路径
        
    Returns:
        dict: {物种名: motif编码}
    """
    motif_data = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # 跳过注释行和空行
                if not line or line.startswith('#'):
                    continue
                
                # 跳过表头
                if line.startswith('基因名') and '二进制编码' in line:
                    continue
                
                # 解析数据行
                if '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        species_name = parts[0].strip()
                        motif_code = parts[1].strip()
                        motif_data[species_name] = motif_code
                        
        print(f"成功解析motif文件: {file_path}")
        print(f"   找到 {len(motif_data)} 个物种的motif编码")
        
        return motif_data
        
    except FileNotFoundError:
        print(f"Motif文件不存在: {file_path}")
        return {}
    except Exception as e:
        print(f"解析motif文件失败: {e}")
        return {}


def merge_data(topology_data, motif_data, verbose=False):
    """
    合并拓扑和motif数据
    
    Args:
        topology_data: 拓扑编码数据
        motif_data: motif编码数据
        verbose: 是否显示详细信息
        
    Returns:
        dict: 合并后的数据
    """
    merged_data = defaultdict(OrderedDict)
    
    for gene_name, species_topo in topology_data.items():
        if verbose:
            print(f"\n处理基因: {gene_name}")
            
        for species_name, topo_code in species_topo.items():
            if species_name in motif_data:
                motif_code = motif_data[species_name]
                merged_data[gene_name][species_name] = {
                    'topology': topo_code,
                    'motif': motif_code
                }
                if verbose:
                    print(f"  {species_name}: {topo_code} | {motif_code}")
            else:
                # 如果没有找到对应的motif编码，用N/A填充
                merged_data[gene_name][species_name] = {
                    'topology': topo_code,
                    'motif': 'N/A'
                }
                if verbose:
                    print(f"  {species_name}: {topo_code} | N/A (motif未找到)")
    
    return merged_data


def write_output(merged_data, output_file, format_style="multi"):
    """
    写入合并结果
    
    Args:
        merged_data: 合并后的数据
        output_file: 输出文件路径
        format_style: 输出格式风格
    """
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            for gene_name, species_data in merged_data.items():
                # 写入基因名
                f.write(f"{gene_name}:\n")
                
                # 写入每个物种的数据
                for species_name, codes in species_data.items():
                    topo_code = codes['topology']
                    motif_code = codes['motif']
                    
                    if format_style == "multi":
                        # 按照multi文件的格式：物种名\t拓扑编码 motif编码
                        f.write(f"{species_name}\t{topo_code} {motif_code}\n")
                    elif format_style == "tab":
                        # 制表符分隔格式：物种名\t拓扑编码\tmotif编码
                        f.write(f"{species_name}\t{topo_code}\t{motif_code}\n")
                    elif format_style == "detailed":
                        # 详细格式
                        f.write(f"{species_name}\t{topo_code}\t{motif_code}\t{len(topo_code)}\t{len(motif_code)}\n")
                
                f.write("\n")  # 基因之间空行分隔
        
        print(f"结果已保存到: {output_file}")
        
    except Exception as e:
        print(f"写入输出文件失败: {e}")


def print_summary(merged_data):
    """
    打印数据摘要
    
    Args:
        merged_data: 合并后的数据
    """
    print("\n" + "="*60)
    print("数据合并摘要")
    print("="*60)
    
    total_genes = len(merged_data)
    total_species = sum(len(species_data) for species_data in merged_data.values())
    
    print(f"总基因数: {total_genes}")
    print(f"总物种-基因组合数: {total_species}")
    
    for gene_name, species_data in merged_data.items():
        print(f"\n基因 {gene_name}:")
        print(f"  物种数: {len(species_data)}")
        
        # 统计编码长度
        topo_lengths = [len(data['topology']) for data in species_data.values() if data['topology'] != 'N/A']
        motif_lengths = [len(data['motif']) for data in species_data.values() if data['motif'] != 'N/A']
        
        if topo_lengths:
            print(f"  拓扑编码长度范围: {min(topo_lengths)}-{max(topo_lengths)}")
        if motif_lengths:
            print(f"  Motif编码长度范围: {min(motif_lengths)}-{max(motif_lengths)}")
        
        # 显示缺失数据
        missing_motif = sum(1 for data in species_data.values() if data['motif'] == 'N/A')
        if missing_motif > 0:
            print(f"    缺失motif编码的物种数: {missing_motif}")


def validate_files(*file_paths):
    """
    验证输入文件是否存在
    
    Args:
        *file_paths: 文件路径列表
        
    Returns:
        bool: 所有文件都存在则返回True
    """
    all_exist = True
    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"文件不存在: {file_path}")
            all_exist = False
        else:
            file_size = os.path.getsize(file_path)
            print(f"文件存在: {file_path} ({file_size} bytes)")
    
    return all_exist


def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(
        description="合并拓扑编码和motif编码文件",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python merge_topology_motif.py -t topology_codes_summary.txt -m gene_motif_binary.txt -o merged_result.txt
  
  # 使用参考格式文件
  python merge_topology_motif.py -t topo.txt -m motif.txt -r multi_new.txt -o result.txt
  
  # 详细模式
  python merge_topology_motif.py -t topo.txt -m motif.txt -o result.txt --verbose
  
  # 使用不同的输出格式
  python merge_topology_motif.py -t topo.txt -m motif.txt -o result.txt --format tab
        """
    )
    
    parser.add_argument('-t', '--topology', 
                       required=True,
                       help='拓扑编码文件路径 (topology_codes_summary.txt)')
    
    parser.add_argument('-m', '--motif',
                       required=True, 
                       help='motif二进制编码文件路径 (gene_motif_binary.txt)')
    
    parser.add_argument('-r', '--reference',
                       help='参考格式文件路径 (multi_new.txt) - 可选')
    
    parser.add_argument('-o', '--output',
                       required=True,
                       help='输出文件路径')
    
    parser.add_argument('--format',
                       choices=['multi', 'tab', 'detailed'],
                       default='multi',
                       help='输出格式 (默认: multi)')
    
    parser.add_argument('--verbose', '-v',
                       action='store_true',
                       help='显示详细处理信息')
    
    args = parser.parse_args()
    
    print("拓扑编码和Motif编码合并工具")
    print("="*50)
    
    # 验证输入文件
    files_to_check = [args.topology, args.motif]
    if args.reference:
        files_to_check.append(args.reference)
        
    if not validate_files(*files_to_check):
        print("\n部分输入文件不存在，请检查文件路径")
        sys.exit(1)
    
    # 解析拓扑文件
    print(f"\n解析拓扑编码文件...")
    topology_data = parse_topology_file(args.topology, args.verbose)
    if not topology_data:
        print("未能解析拓扑编码文件")
        sys.exit(1)
    
    # 解析motif文件
    print(f"\n解析motif编码文件...")
    motif_data = parse_motif_file(args.motif)
    if not motif_data:
        print("未能解析motif编码文件")
        sys.exit(1)
    
    # 合并数据
    print(f"\n合并数据...")
    merged_data = merge_data(topology_data, motif_data, args.verbose)
    
    if not merged_data:
        print("未能合并任何数据")
        sys.exit(1)
    
    # 写入输出文件
    print(f"\n写入输出文件...")
    write_output(merged_data, args.output, args.format)
    
    # 显示摘要
    if args.verbose:
        print_summary(merged_data)
    
    # 提供使用建议
    print(f"\n输出格式说明:")
    if args.format == 'multi':
        print("  - 基因名:")
        print("  - 物种名<tab>拓扑编码<space>motif编码")
    elif args.format == 'tab':
        print("  - 基因名:")
        print("  - 物种名<tab>拓扑编码<tab>motif编码") 
    elif args.format == 'detailed':
        print("  - 基因名:")
        print("  - 物种名<tab>拓扑编码<tab>motif编码<tab>拓扑长度<tab>motif长度")
    
    print("\n处理完成!")


if __name__ == "__main__":
    main()