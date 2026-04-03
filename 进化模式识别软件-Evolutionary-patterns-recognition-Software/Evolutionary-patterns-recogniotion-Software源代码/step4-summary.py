#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
进化模式识别引擎 - 完整版
处理规范化拓扑编码和motif数量信息
支持外部参数输入
"""

import numpy as np
import re
import argparse
import sys

# 设置程序输出的编码方式为utf-8
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def parse_multi_txt(file_path):
    """
    解析multi_new.txt文件，返回结构化数据
    """
    data = {}
    current_gene = None
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                # 检查是否是基因标题行
                if line.endswith(':'):
                    current_gene = line[:-1]  # 移除冒号
                    data[current_gene] = {}
                elif current_gene is not None:
                    # 使用正则表达式分割，处理制表符和空格混合的情况
                    parts = re.split(r'\s+', line, maxsplit=1)
                    if len(parts) < 2:
                        continue
                        
                    species = parts[0]
                    # 合并所有编码部分(移除空格)
                    code = ''.join(parts[1:]).replace(" ", "")
                    data[current_gene][species] = code
                    
        return data
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        return {}
    except Exception as e:
        print(f"解析文件时出错: {e}")
        return {}

def normalize_codes(species_codes):
    """
    标准化编码：找到最大长度，用占位符X填充较短编码的前面
    """
    # 找到最大编码长度
    max_len = max(len(code) for code in species_codes.values())
    
    normalized = {}
    for species, code in species_codes.items():
        # 在前面添加X占位符
        padding = "X" * (max_len - len(code))
        normalized[species] = padding + code
    
    return normalized

class TrieNode:
    """前缀树节点类"""
    def __init__(self):
        self.children = {}  # 子节点字典
        self.is_end = False  # 是否代表一个编码的结束
        self.species = []  # 存储到达此节点的物种名称

class Trie:
    """前缀树类"""
    def __init__(self):
        self.root = TrieNode()
    
    def insert(self, species, code):
        """
        向前缀树中插入一个物种及其编码
        """
        node = self.root
        for char in code:
            # 跳过占位符X
            if char == 'X':
                continue
                
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
            # 在路径上的每个节点都记录物种
            node.species.append(species)
        node.is_end = True
    
    def find_clusters(self, min_cluster_size=2):
        """
        查找所有包含至少min_cluster_size个物种的簇
        """
        clusters = []
        
        def dfs(node, current_prefix):
            # 如果当前节点包含足够多的物种，记录为一个簇
            if len(node.species) >= min_cluster_size:
                clusters.append({
                    'prefix': current_prefix,
                    'species': node.species.copy(),  # 使用副本避免后续修改
                    'size': len(node.species)
                })
            
            # 递归遍历所有子节点
            for char, child_node in node.children.items():
                dfs(child_node, current_prefix + char)
        
        # 从根节点开始深度优先搜索
        dfs(self.root, "")
        
        # 按前缀长度排序(最长的在前)
        clusters.sort(key=lambda x: len(x['prefix']), reverse=True)
        return clusters

def weighted_hamming_distance(s1, s2):
    """
    计算两个字符串的加权汉明距离
    处理不同长度和motif数量信息(数字2等)
    """
    # 处理不同长度的情况，取较短的长度
    min_len = min(len(s1), len(s2))
    distance = 0
    
    for i in range(min_len):
        c1, c2 = s1[i], s2[i]
        
        # 跳过占位符
        if c1 == 'X' or c2 == 'X':
            continue
            
        # 处理motif数量差异
        if c1 != c2:
            # 如果都是数字，计算数值差异
            if c1.isdigit() and c2.isdigit():
                num1, num2 = int(c1), int(c2)
                distance += abs(num1 - num2)
            else:
                # 其他情况，差异为1
                distance += 1
    
    # 加上长度差异
    distance += abs(len(s1) - len(s2))
    
    return distance

def calculate_cluster_stats(cluster_species, species_codes):
    """
    计算簇的统计信息
    """
    if len(cluster_species) < 2:
        return 0, 0, []
    
    distances = []
    codes = [species_codes[sp] for sp in cluster_species]
    
    # 计算所有两两组合的距离
    for i in range(len(codes)):
        for j in range(i+1, len(codes)):
            dist = weighted_hamming_distance(codes[i], codes[j])
            distances.append(dist)
    
    if distances:
        mu = np.mean(distances)
        sigma = np.std(distances)
        return mu, sigma, distances
    else:
        return 0, 0, []

def detect_anomalies(cluster, species_codes, z_threshold_high=3, z_threshold_low=-2):
    """
    检测簇中的异常个体
    """
    cluster_species = cluster['species']
    mu, sigma, _ = calculate_cluster_stats(cluster_species, species_codes)
    anomalies = {'high': [], 'low': []}
    
    if sigma == 0 or len(cluster_species) < 2:
        return anomalies
    
    for species in cluster_species:
        # 计算该物种到簇内所有其他物种的平均距离
        dists = []
        for other in cluster_species:
            if other != species:
                dist = weighted_hamming_distance(species_codes[species], species_codes[other])
                dists.append(dist)
        
        if not dists:
            continue
            
        mean_dist = np.mean(dists)
        z_score = (mean_dist - mu) / sigma if sigma != 0 else 0
        
        if z_score > z_threshold_high:
            anomalies['high'].append({
                'species': species,
                'z_score': z_score,
                'mean_distance': mean_dist
            })
        elif z_score < z_threshold_low:
            anomalies['low'].append({
                'species': species,
                'z_score': z_score,
                'mean_distance': mean_dist
            })
    
    return anomalies

def analyze_gene(gene_name, species_codes, min_cluster_size=2, z_threshold_high=3, z_threshold_low=-2):
    """
    分析单个基因的数据
    """
    print(f"\n=== 分析基因: {gene_name} ===")
    print(f"物种数量: {len(species_codes)}")
    
    # 标准化编码长度
    normalized_codes = normalize_codes(species_codes)
    
    # 打印所有编码以进行调试
    print("标准化后的编码样本:")
    for i, (species, code) in enumerate(list(normalized_codes.items())[:3]):
        print(f"  {species}: {code} (原始: {species_codes[species]})")
    if len(normalized_codes) > 3:
        print("  ...")
    
    # 构建前缀树
    trie = Trie()
    for species, code in normalized_codes.items():
        trie.insert(species, code)
    
    # 查找簇
    clusters = trie.find_clusters(min_cluster_size)
    print(f"找到 {len(clusters)} 个簇 (最小大小: {min_cluster_size})")
    
    results = {
        'gene': gene_name,
        'clusters': [],
        'anomalies': []
    }
    
    # 分析每个簇
    for i, cluster in enumerate(clusters):
        print(f"\n簇 {i+1}:")
        print(f"  公共前缀: '{cluster['prefix']}' (长度: {len(cluster['prefix'])})")
        print(f"  包含物种: {', '.join(cluster['species'])}")
        print(f"  簇大小: {cluster['size']}")
        
        # 计算簇统计量
        mu, sigma, distances = calculate_cluster_stats(cluster['species'], normalized_codes)
        print(f"  平均加权汉明距离: {mu:.2f} ± {sigma:.2f}")
        
        # 检测异常
        anomalies = detect_anomalies(cluster, normalized_codes, z_threshold_high, z_threshold_low)
        
        if anomalies['high']:
            print("  异常个体 (Z > 3, 可能特殊进化事件):")
            for anomaly in anomalies['high']:
                print(f"    * {anomaly['species']}: Z值={anomaly['z_score']:.2f}, 平均距离={anomaly['mean_distance']:.2f}")
        
        if anomalies['low']:
            print("  异常个体 (Z < -2, 可能进化停滞):")
            for anomaly in anomalies['low']:
                print(f"    * {anomaly['species']}: Z值={anomaly['z_score']:.2f}, 平均距离={anomaly['mean_distance']:.2f}")
        
        # 保存结果
        cluster_info = {
            'id': i+1,
            'prefix': cluster['prefix'],
            'species': cluster['species'],
            'size': cluster['size'],
            'mean_distance': mu,
            'std_distance': sigma,
            'anomalies': anomalies
        }
        results['clusters'].append(cluster_info)
        
        # 收集所有异常
        results['anomalies'].extend(anomalies['high'])
        results['anomalies'].extend(anomalies['low'])
    
    return results

def detect_hgt(species, all_codes, topology_codes, hgt_threshold=5):
    """
    检测水平基因转移(HGT)
    基于拓扑编码和特征编码的矛盾
    """
    candidate = species
    candidate_code = all_codes[candidate]
    
    # 1. 找到拓扑编码最相似的物种(理论近邻)
    min_topology_dist = float('inf')
    theoretical_neighbor = None
    
    for other, topo_code in topology_codes.items():
        if other != candidate:
            dist = weighted_hamming_distance(topo_code, topology_codes[candidate])
            if dist < min_topology_dist:
                min_topology_dist = dist
                theoretical_neighbor = other
    
    if theoretical_neighbor is None:
        return False
    
    # 2. 计算与理论近邻的特征差异
    theoretical_neighbor_code = all_codes[theoretical_neighbor]
    theoretical_diff = weighted_hamming_distance(candidate_code, theoretical_neighbor_code)
    
    # 3. 找到实际在特征空间中的最近邻
    min_feature_dist = float('inf')
    actual_neighbor = None
    
    for other, code in all_codes.items():
        if other != candidate:
            dist = weighted_hamming_distance(candidate_code, code)
            if dist < min_feature_dist:
                min_feature_dist = dist
                actual_neighbor = other
    
    # 4. 应用逻辑判断
    if theoretical_neighbor != actual_neighbor and (min_feature_dist + hgt_threshold) < theoretical_diff:
        print(f"HGT候选: {candidate}")
        print(f"  拓扑近邻: {theoretical_neighbor} (距离: {min_topology_dist})")
        print(f"  特征近邻: {actual_neighbor} (距离: {min_feature_dist})")
        print(f"  与拓扑近邻的特征差异: {theoretical_diff}")
        return True
    
    return False

def detect_accelerated_evolution(cluster, species_codes):
    """
    检测加速进化
    基于特征编码长度异常
    """
    accelerated_species = []
    
    # 计算簇内平均编码长度(忽略占位符)
    lengths = [len(code.replace('X', '')) for code in species_codes.values()]
    avg_length = np.mean(lengths)
    
    for species, code in species_codes.items():
        # 计算有效编码长度(忽略占位符)
        effective_length = len(code.replace('X', ''))
        
        # 检查是否显著长于平均值
        if effective_length > 1.5 * avg_length:
            accelerated_species.append({
                'species': species,
                'length': effective_length,
                'avg_length': avg_length
            })
    
    return accelerated_species

def detect_convergent_evolution(species_codes, topology_codes, similarity_threshold=0.7, topology_diff_threshold=3):
    """
    检测趋同进化
    基于远缘物种特征相似
    """
    convergent_pairs = []
    species_list = list(species_codes.keys())
    
    for i in range(len(species_list)):
        for j in range(i+1, len(species_list)):
            sp1, sp2 = species_list[i], species_list[j]
            
            # 计算拓扑距离
            topology_dist = weighted_hamming_distance(topology_codes[sp1], topology_codes[sp2])
            
            # 如果拓扑距离大，检查特征相似度
            if topology_dist > topology_diff_threshold:
                # 计算特征相似度
                feature_sim = 1 - (weighted_hamming_distance(species_codes[sp1], species_codes[sp2]) / 
                                  max(len(species_codes[sp1].replace('X', '')), 
                                      len(species_codes[sp2].replace('X', ''))))
                
                if feature_sim > similarity_threshold:
                    convergent_pairs.append({
                        'species1': sp1,
                        'species2': sp2,
                        'topology_distance': topology_dist,
                        'feature_similarity': feature_sim
                    })
    
    return convergent_pairs

def main():
    """主函数"""
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='进化模式识别引擎')
    parser.add_argument('-i', '--input', type=str, default='/mnt/f/BGI_working/Project_基因组变异和BCR,TCR/Sequence_Alignment/9_species/TRBC/result.txt',
                        help='输入文件路径')
    parser.add_argument('-m', '--min_cluster_size', type=int, nargs='+', default=[2, 3],
                        help='最小簇大小，可以是一个或多个值')
    parser.add_argument('-zh', '--z_threshold_high', type=float, default=3.0,
                        help='高Z值阈值，用于检测特殊进化事件')
    parser.add_argument('-zl', '--z_threshold_low', type=float, default=-2.0,
                        help='低Z值阈值，用于检测进化停滞')
    parser.add_argument('-ht', '--hgt_threshold', type=float, default=5.0,
                        help='水平基因转移检测阈值')
    parser.add_argument('-cs', '--converge_similarity', type=float, default=0.7,
                        help='趋同进化相似度阈值')
    parser.add_argument('-ct', '--converge_topology', type=float, default=3.0,
                        help='趋同进化拓扑差异阈值')
    
    # 解析参数
    args = parser.parse_args()
    
    # 文件路径
    file_path = args.input
    
    print("开始解析数据文件...")
    data = parse_multi_txt(file_path)
    
    if not data:
        print("未能解析出任何数据，请检查文件路径和格式")
        return
    
    print(f"成功解析 {len(data)} 个基因的数据")
    
    # 分析每个基因
    all_results = {}
    for gene_name, species_codes in data.items():
        # 尝试不同的最小簇大小
        for min_size in args.min_cluster_size:
            print(f"\n尝试最小簇大小: {min_size}")
            results = analyze_gene(gene_name, species_codes, min_cluster_size=min_size, 
                                  z_threshold_high=args.z_threshold_high, 
                                  z_threshold_low=args.z_threshold_low)
            if results['clusters']:  # 如果找到了簇
                all_results[gene_name] = results
                break
    
    # 输出总结
    print("\n" + "="*50)
    print("分析总结:")
    print("="*50)
    
    for gene_name, results in all_results.items():
        total_anomalies = len(results['anomalies'])
        print(f"{gene_name}: {len(results['clusters'])} 个簇, {total_anomalies} 个异常个体")
        
        if total_anomalies > 0:
            print("  异常个体:")
            for anomaly in results['anomalies']:
                anomaly_type = "特殊事件" if anomaly['z_score'] > args.z_threshold_high else "进化停滞"
                print(f"    * {anomaly['species']} ({anomaly_type}): Z={anomaly['z_score']:.2f}")
        else:
            print("  无异常个体检测到")
    
    # 高级进化模式检测
    print("\n" + "="*50)
    print("高级进化模式检测:")
    print("="*50)
    
    for gene_name, species_codes in data.items():
        print(f"\n--- 检测基因 {gene_name} 的高级进化模式 ---")
        
        # 标准化编码
        normalized_codes = normalize_codes(species_codes)
        
        # 提取拓扑编码部分 (假设前几位是拓扑编码)
        # 这里需要根据实际情况调整，或者从单独的文件中读取拓扑编码
        topology_codes = {}
        for species, code in normalized_codes.items():
            # 提取拓扑编码部分 (这里假设前5位是拓扑编码，需要根据实际情况调整)
            topology_part = code.replace('X', '')[:5]  # 调整这个数字
            topology_codes[species] = topology_part
        
        # 检测水平基因转移
        print("水平基因转移(HGT)检测:")
        hgt_candidates = []
        for species in species_codes.keys():
            if detect_hgt(species, normalized_codes, topology_codes, args.hgt_threshold):
                hgt_candidates.append(species)
        
        if not hgt_candidates:
            print("  未发现HGT候选")
        
        # 检测加速进化
        print("\n加速进化检测:")
        accelerated = detect_accelerated_evolution(results['clusters'][0] if results['clusters'] else {}, normalized_codes)
        if accelerated:
            for acc in accelerated:
                print(f"  {acc['species']}: 编码长度={acc['length']:.1f}, 平均长度={acc['avg_length']:.1f}")
        else:
            print("  未发现加速进化候选")
        
        # 检测趋同进化
        print("\n趋同进化检测:")
        convergent = detect_convergent_evolution(normalized_codes, topology_codes, 
                                                args.converge_similarity, args.converge_topology)
        if convergent:
            for conv in convergent:
                print(f"  {conv['species1']} 和 {conv['species2']}: "
                      f"拓扑距离={conv['topology_distance']}, 特征相似度={conv['feature_similarity']:.2f}")
        else:
            print("  未发现趋同进化候选")
    
    print("\n分析完成!")

if __name__ == "__main__":
    main()