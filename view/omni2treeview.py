import csv
import argparse
from Bio import Phylo
import csv
import json
import argparse
from collections import defaultdict
from typing import Dict, Any
from typing import Tuple
import logging



def parse_newick_to_csv(nwk_file: str):

    out_dict = {
        "header": ["parent", "node", "branch.length","label","confidence"],
        "rows": [],
        "rows_final": []
    }

    unique_names = set()

    tree = Phylo.read(nwk_file, 'newick')
    
    internal_counter = 0
    node_names = {}  
    
    for clade in tree.find_clades():
        if clade.name:
            node_names[id(clade)] = clade.name
        else:
            internal_counter += 1
            node_names[id(clade)] = f'Internal-{internal_counter}'



    if tree.root.name is None:
        internal_counter += 1
        tree.root.name = f'Internal-{internal_counter}'
        node_names[id(tree.root)] = tree.root.name
    
    def clean_leaf_name(name: str) -> str:
        """ leaf name has strucure like prefix_name_suffix, we only keep name part """
        nn = name.split("_")
        # print(nn)
        if len(nn) > 1:
            return nn[1]
        else:
            return name
        
    def traverse(clade, parent_name):
        node_name = node_names[id(clade)]
        node_name = clean_leaf_name(node_name)
        branch_length = clade.branch_length if clade.branch_length is not None else 0.0
        out_dict["rows"].append([parent_name, node_name, branch_length, node_name, "NA"])
        unique_names.add(parent_name)
        unique_names.add(node_name)
        
        for child in clade.clades:
            traverse(child, node_name)
    root_name = node_names[id(tree.root)]
    traverse(tree.root, root_name)

    uniq_name_2_idx = {}
    idx = 0
    for uname in unique_names:
        uniq_name_2_idx[uname] = idx
        idx += 1

    # process each row from 'rows' to 'rows_final', assign unique index to each node,
    # put the node_name to the last column
    for row in out_dict["rows"]:
        parent_name = row[0]
        node_name = row[1]
        branch_length = row[2]
        parent_idx = uniq_name_2_idx[parent_name]
        node_idx = uniq_name_2_idx[node_name]
        out_dict["rows_final"].append([parent_idx, node_idx, branch_length, node_name, "NA"])
    

    return out_dict

def load_meta(meta_file: str):
    meta_dict = {
        "header": [],
        "col_type": [],
        "rows": []
    }
    with open(meta_file, 'r') as f:
        reader = csv.DictReader(f)
        uniq_ids = set()
        meta_dict["header"] = reader.fieldnames
        # first row to determine column types
        first_row = next(reader)
        for key, value in first_row.items():
            if value.lower() in ['character', 'date','numeric','integer']:
                meta_dict["col_type"].append(value.lower())
            else:
                logging.error(f"Unknown column type: {value} for column {key}. Please check the column type in the second row, eligible values are [character, date, numeric, integer]. Exiting...")
                exit(1)

        # process next rows
        for row in reader:
            # get the first column as sample_id
            sample_id = row[meta_dict["header"][0]]
            # sample_id_list = sample_id.split('_')
            # sample_id = sample_id_list[1]
            if sample_id in uniq_ids:
                logging.error(f"Duplicate sample_id {sample_id} found in metadata. exiting...")
                exit(1)
            uniq_ids.add(sample_id)

            # fill in missing columns with NA
            for key in meta_dict["header"]:
                if key not in row or row[key] == '':
                    row[key] = 'NA'

            meta_dict[sample_id] = row

    # validate the values based on col_type
    for sample_id, row in meta_dict.items():
        if sample_id in ['header', 'col_type', 'rows']:
            continue
        for idx, key in enumerate(meta_dict["header"]):
            value = row[key]
            if value == 'NA':
                continue
            col_type = meta_dict["col_type"][idx]
            if col_type == 'numeric' or col_type == 'integer':
                try:
                    if col_type == 'numeric':
                        float(value)
                    else:
                        int(value)
                except ValueError:
                    logging.error(f"Value '{value}' for column '{key}' in sample_id '{sample_id}' is not of type {col_type}. exiting...")
                    exit(1)
            elif col_type == 'date':
                # simple check for date format (could be improved)
                if not value.count('/') == 2 and not value.count('-') == 2 and not value.count('-') == 1:

                    if not value.isdigit():

                        logging.error(f"Value '{value}' for column '{key}' in sample_id '{sample_id}' is not a valid Date format. exiting...")
                        exit(1)
            # character type does not need validation
    # print(meta_dict)
    return meta_dict

def load_code(code_file: str):


    code_dict = {}
    
    if code_file is None:
        return code_dict
    
    internal_acc_set = set()
    with open(code_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                logging.error(f"Invalid line in code file: {line}. exiting...")
                exit(1)
            acc, internal_acc = parts
            if internal_acc in internal_acc_set:
                logging.error(f"Duplicate internal accession {internal_acc} found in code file. exiting...")
                exit(1)
            internal_acc_set.add(internal_acc)
            code_dict[internal_acc] = acc
    return code_dict

def parse_csv_to_tree(input_file: str,tree_name: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    nodes = {}
    internal_node_count = 0
    leaf_node_count = 0
    children_map = defaultdict(list)
    # meta_keys = {}

    common_header = ["parent", "node", "branch.length", "label","confidence"]
    common_header = [ h.capitalize() for h in common_header]

    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        types = next(reader)
        # print("Detected column types:")
        # for k, v in types.items():
        #     print(f"{k}: {v}")
        if not all(t in ["numeric", "integer", "date","character"] for t in types.values()):
            raise ValueError("First row of CSV must contain type definitions like 'numeric', 'integer', or 'string'.")

        meta_keys = {k:{
            "type":v,
            "std_type":None,
            "name":k,
            "NA_count":0,
            "min":float("inf"),
            "max": float("-inf"),
            "categories":[]
        }  for k, v in types.items() if k not in common_header}

        # print(meta_keys)

        for row in reader:
            node_id = row["Node"]
            parent_id = row["Parent"]

            meta = {
                k: row[k]
                for k in row
                if k not in common_header
            }
            # print(row)
            # print(meta)
        

            for k, v in meta.items():
                if meta_keys[k]["type"] in ["numeric", "integer"]: # type: ignore
                    if v == "NA":
                        meta_keys[k]["NA_count"] += 1
                    else:
                        try:
                            meta_keys[k]["min"] = min(meta_keys[k]["min"], float(v))
                            meta_keys[k]["max"] = max(meta_keys[k]["max"], float(v))
                        except ValueError:
                            logging.warning(f"Non-numeric value '{v}' found in numeric column '{k}'. Treating as NA.")
                            meta_keys[k]["NA_count"] += 1
                    if meta_keys[k]["std_type"] is None:
                        meta_keys[k]["std_type"] = "cont"
                else:
                    if v == "NA":
                        meta_keys[k]["NA_count"] += 1
                    else:
                        if v not in meta_keys[k]["categories"]:
                            meta_keys[k]["categories"].append(v)
                    if meta_keys[k]["std_type"] is None:
                        meta_keys[k]["std_type"] = "cate"

            brlen = 0.0 if row["Branch.length"] == "NA" else float(row["Branch.length"])

            nodes[node_id] = {
                "id": node_id,
                "name": row["Label"],
                "branch_length": brlen,
                "meta": meta,
                "confidence": row["Confidence"],
                "children": []
            }

            if parent_id != node_id:
                children_map[parent_id].append(node_id)
            else:
                root_id = node_id
        
            # meta_keys = list(meta.keys())
        # print("Meta keys:", json.dumps(meta_keys, indent=2))
        


    def build_tree(node_id: str) -> Dict[str, Any]:
        node = nodes[node_id]
        node["children"] = [build_tree(child_id) for child_id in children_map.get(node_id, [])]
        return node
    
    if 'root_id' not in locals():
        raise ValueError("Root node not found. Ensure there is a row where 'parent' == 'node' to define the root.")

    tree = build_tree(root_id)


    for node_id, node in nodes.items():
        if len(node['children']) == 0:
            leaf_node_count += 1
        else:
            internal_node_count += 1

    tree_meta = {
        "name": tree_name,
        "node_count": len(nodes),
        "leaf_count": leaf_node_count,
        "internal_count": internal_node_count,
        "meta_keys": meta_keys
    }

    return  (tree, tree_meta)

def estimate_tree_dimensions(tree: Dict[str, Any]) -> Dict[str, int]:
    def dfs(node, depth):
        nonlocal max_depth, leaf_count
        max_depth = max(max_depth, depth)
        if not node.get("children"):  # 叶子节点
            leaf_count += 1
        for child in node.get("children", []):
            dfs(child, depth + 1)

    max_depth = 0
    leaf_count = 0
    dfs(tree, 1)

    return {
        "estimated_width": max_depth,
        "estimated_height": leaf_count
    }

def main(input_file: str, out_prefix: str, html_template: str, tree_name: str):
    tree,tree_meta = parse_csv_to_tree(input_file, tree_name)


    # clean the inf and -inf values

    for k, v in tree_meta["meta_keys"].items():
        if v["min"] == float("inf"):
            v["min"] = None
        if v["max"] == float("-inf"):
            v["max"] = None

    # estimate the tree dimensions
    dimensions = estimate_tree_dimensions(tree)
    tree_meta["estimated_width"] = dimensions["estimated_width"]
    tree_meta["estimated_height"] = dimensions["estimated_height"]


    out_html = out_prefix + ".html"
    tree_json_output = out_prefix + ".tree.json"
    tree_meta_json_output = out_prefix + ".tree_meta.json"
    with open(tree_json_output, "w") as f:
        json.dump(tree, f, indent=2)
    with open(tree_meta_json_output, "w") as f2:
        json.dump(tree_meta, f2, indent=2)
    tree_str = json.dumps(tree)
    with open(html_template, "r") as f:
        template = f.read()
    with open(out_html, "w") as f2:

        template = template.replace("[/*inject_meta*/]", json.dumps(tree_meta))

        f2.write(template.replace("[/*inject_data*/]", tree_str))




if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(description='Convert Newick tree to CSV format.')
    parser.add_argument('-n', '--nwk_file', type=str, help='Input Newick file')
    parser.add_argument('-m', '--meta', type=str, help='Metadata file')
    parser.add_argument('-c', '--code', type=str, help='conversion code file',default=None)
    parser.add_argument('-o', '--output', type=str, help='Output prefix, without extension')
    parser.add_argument("-t","--template",required=True, help="Path to input HTML template")
    parser.add_argument("-l","--name",required=True, help="Name of the tree",default="omnitreeview")

    args = parser.parse_args()

    # print arguments for logging
    logging.info("> Arguments:")
    for arg in vars(args):
        logging.info(f"  {arg}: {getattr(args, arg)}")

    # handle the output format, if it has any extension, remove it
    output_base = args.output
    
    if output_base == '.' or output_base.endswith('./') or output_base.endswith('.\\'):
        logging.error("Error: Output prefix cannot end with a dot or a path separator. exiting...")
        exit(1)
    
    if '.' in output_base:
        logging.warning(f"Output prefix '{output_base}' contains an extension. Removing the extension for output base name.")
        output_base = '.'.join(output_base.split('.')[:-1])
    
    
    tree_csv = parse_newick_to_csv(args.nwk_file)

    meta_dict = load_meta(args.meta)

    code_dict = load_code(args.code)


    final_header = tree_csv["header"] + meta_dict["header"][1:]  # skip sample_id column in meta
    final_header = [x.capitalize() for x in final_header]
    final_col_type = ['integer','integer','numeric','character','numeric'] + meta_dict["col_type"][1:]
    logging.info(f"> Final header and column types ({len(final_header)}):")
    for idx, col in enumerate(final_header):
        logging.info(f"  {col} --> {final_col_type[idx]}")
    
    out_fh = open(output_base+".meta.csv", 'w', newline='')
    writer = csv.writer(out_fh, lineterminator='\n')
    writer.writerow(final_header)   
    writer.writerow(final_col_type)

    # print(meta_dict["header"])

    logging.info(f"> Generating merged CSV with metadata to {output_base}.meta.csv")
    all_meta_sample_ids = set(meta_dict.keys()) - set(["header", "col_type", "rows"])
    for row in tree_csv["rows_final"]:
        node_name = row[3]
        # first try check meta dict, then code dict
        if node_name in meta_dict:
            meta_row = meta_dict[node_name]
            all_meta_sample_ids.discard(node_name)
        elif node_name in code_dict and code_dict[node_name] in meta_dict:
            meta_row = meta_dict[code_dict[node_name]]
            all_meta_sample_ids.discard(code_dict[node_name])
        else:
            # print warnings
            if not node_name.startswith("Internal-"):
                logging.warning(f"No metadata found for node '{node_name}'. Filling with NA. If this is not intentional, please check your metadata file and code file for consistency with the tree node names(second field after splitted by '_').")
            # fill with NA
            meta_row = {key: 'NA' for key in meta_dict["header"]}
        # remove sample_id column
        meta_values = [meta_row[key] for key in meta_dict["header"][1:]]
        # print(row)
        final_row = row + meta_values
        writer.writerow(final_row)
    if len(all_meta_sample_ids) > 0:
        logging.warning(f"There are {len(all_meta_sample_ids)} sample IDs in metadata that are not found in the tree nodes. These samples will be ignored. Sample IDs: {', '.join(all_meta_sample_ids)}, please check for consistency between your metadata file and tree node names.")
    out_fh.close()

    logging.info("> Generating tree and HTML output...")
    main(output_base+".meta.csv", output_base, args.template, args.name)
    logging.info("> Done.")
    logging.info("> Outputs generated:")
    logging.info(f" {output_base}.meta.csv")
    logging.info(f" {output_base}.tree.json")
    logging.info(f" {output_base}.tree_meta.json")
    logging.info(f" {output_base}.html")