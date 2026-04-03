import os
import argparse
from ete3 import Tree

def assign_topology_code(tree):
    """
    为树的每个叶子节点分配拓扑编码
    """
    code_map = {}
    
    def traverse(node, path=""):
        if node.is_leaf():
            # 保存叶子节点的编码
            code_map[node.name] = path
        else:
            # 获取所有子节点
            children = node.children
            
            # 确保是二叉结构
            if len(children) != 2:
                print(f"[警告] 在节点 {node.name} 发现非二叉结构 ({len(children)} 个子节点)")
                # 尝试处理多叉树：按顺序分配0,1,2...
                for idx, child in enumerate(children):
                    traverse(child, path + str(idx))
            else:
                # 二叉结构：第一个孩子为左（0），第二个为右（1）
                traverse(children[0], path + "0")  # 左分支
                traverse(children[1], path + "1")  # 右分支

    traverse(tree, "")
    return code_map

def read_newick_and_encode(file_path):
    """
    读取 Newick 文件，构建树，并生成拓扑编码
    """
    try:
        with open(file_path, 'r', encoding="utf-8") as f:
            newick_str = f.read().strip()

        # 解析 Newick 字符串（忽略分支长度和置信度标签）
        tree = Tree(newick_str, format=1)  # format=1 支持带标签和长度的格式
        return assign_topology_code(tree)
    except Exception as e:
        print(f"[ERROR] 无法解析文件 {file_path}: {e}")
        return {}

def process_files(input_paths, output_dir):
    """
    处理输入的 Newick 文件并生成拓扑编码
    """
    all_results = {}
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    for gene, path in input_paths.items():
        if not os.path.exists(path):
            print(f"[警告] 文件不存在: {path}")
            continue

        print(f"处理基因: {gene}")
        codes = read_newick_and_encode(path)
        all_results[gene] = codes
        
        # 保存每个基因的单独结果
        output_file = os.path.join(output_dir, f"{gene}_topology_codes.txt")
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(f"Gene: {gene}\n")
            f.write("物种名称\t拓扑编码\n")
            for leaf, code in codes.items():
                f.write(f"{leaf}\t{code}\n")
        print(f"  结果已保存至: {output_file}")
        print(f"  找到 {len(codes)} 个叶子节点")
        print("-" * 50)

    # 保存所有结果到汇总文件
    if all_results:
        summary_file = os.path.join(output_dir, "topology_codes_summary.txt")
        with open(summary_file, "w", encoding="utf-8") as out_f:
            for gene, codes in all_results.items():
                out_f.write(f"Gene: {gene}\n")
                out_f.write("物种名称\t拓扑编码\n")
                for leaf, code in codes.items():
                    out_f.write(f"{leaf}\t{code}\n")
                out_f.write("\n" + "="*50 + "\n\n")
        print(f"所有结果已保存至汇总文件: {summary_file}")
    else:
        print("没有找到任何有效的Newick文件或无法解析任何文件")

def main():
    # 安装必要的依赖（如果尚未安装）
    try:
        from ete3 import Tree
    except ImportError:
        print("安装必要的Python库...")
        import subprocess
        subprocess.run(["pip", "install", "ete3", "numpy", "six"], check=True)
        print("安装完成，请重新运行脚本")
        exit(0)
    
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='为系统发育树生成拓扑编码')
    
    # 添加参数
    parser.add_argument('-i', '--input', nargs='+', required=True,
                       help='输入文件路径，可以是一个或多个Newick格式文件')
    parser.add_argument('-o', '--output', default='./output',
                       help='输出目录路径，默认为当前目录下的output文件夹')
    parser.add_argument('-n', '--names', nargs='+',
                       help='为每个输入文件指定基因名称，数量应与输入文件相同')
    
    # 解析参数
    args = parser.parse_args()
    
    # 处理输入文件
    input_files = args.input
    
    # 处理基因名称
    if args.names:
        if len(args.names) != len(input_files):
            print("错误: 提供的基因名称数量与输入文件数量不匹配")
            exit(1)
        gene_names = args.names
    else:
        # 如果没有提供基因名称，使用文件名（不含扩展名）作为基因名
        gene_names = [os.path.splitext(os.path.basename(f))[0] for f in input_files]
    
    # 创建文件路径字典
    file_paths = {name: path for name, path in zip(gene_names, input_files)}
    
    print("开始处理 Newick 文件...\n")
    print(f"输入文件: {', '.join(input_files)}")
    print(f"输出目录: {args.output}")
    print(f"基因名称: {', '.join(gene_names)}")
    print("-" * 50)
    
    # 处理文件
    process_files(file_paths, args.output)

if __name__ == "__main__":
    main()