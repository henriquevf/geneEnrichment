import argparse
import pandas as pd
from gprofiler import GProfiler

def read_gene_list(file_path):
    with open(file_path, "r") as file:
        genes = [line.strip() for line in file.readlines()]
    return genes

def run_gprofiler(genes, output_prefix):
    gp = GProfiler(return_dataframe=True)

    enrichment_results = gp.profile(organism='hsapiens', query=genes)

    # Save the raw results to a CSV file
    enrichment_results.to_csv(f"{output_prefix}_raw_enrichment_results.csv", index=False)

    return enrichment_results

def parse_and_save_results(enrichment_results, output_prefix):
    filtered_results = enrichment_results[
        (enrichment_results["p_value"] <= 0.05) &
        (enrichment_results["source"].isin(["GO:BP", "GO:MF", "GO:CC", "KEGG", "HP"]))
    ]
    
    sorted_results = filtered_results.sort_values("p_value")

    # Save all significant categories to a new CSV file
    sorted_results.to_csv(f"{output_prefix}_significant_enrichment_categories.csv", index=False)

    # Save all significant categories to a new Excel file
    sorted_results.to_excel(f"{output_prefix}_significant_enrichment_categories.xlsx", index=False)

    return sorted_results

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Run g:Profiler and parse enrichment analysis results")
parser.add_argument("gene_list_file", help="Input file containing the list of genes")
parser.add_argument("output_prefix", help="Output file prefix for significant enrichment categories")

args = parser.parse_args()

# Read the gene list from the input file
genes = read_gene_list(args.gene_list_file)

# Run g:Profiler and save the raw results
enrichment_results = run_gprofiler(genes, args.output_prefix)

# Parse the results and save significant categories
significant_categories = parse_and_save_results(enrichment_results, args.output_prefix)

# Display significant categories
print(significant_categories)

