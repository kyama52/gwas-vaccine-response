"""
# Sub routine script for OrdinalGWAS
# Author: K Yamazaki
# Created: Jan 31 2025

# Change Log:
# Jan 31 2025   Initial release
# Feb 04 2025
#   - Change test method from score test to likelihood ratio test (LRT)
"""

using OrdinalGWAS
using CSV, DataFrames

phenotype = ENV["PHENOTYPE"]
cov_list = ENV["COV_LIST"]
ord_file = ENV["ORD_FILE"]
bim_name = ENV["BIM_NAME"]
fam_file = bim_name * ".fam"
out_prefix = ENV["OUT_PREFIX"]

# Determine whether the delimiter is a space or a tab
first_line = readline(fam_file)
delim = occursin("\t", first_line) ? "\t" : " "

# Load covariate and phenotype file
cov_df = CSV.read(ord_file, DataFrame)
fam_df = CSV.read(fam_file, DataFrame, header=false, delim=delim)
plink_sample_ids = fam_df[:, 2]

# Extract index in plink file from ord_file
sample_ids = cov_df[:, "iid"]
sample_indices = findall(in(sample_ids), plink_sample_ids)

formula_expr = @eval @formula($(Symbol(phenotype)) ~ $(Meta.parse(cov_list)))

ordinalgwas(
    formula_expr,
    ord_file,
    bim_name,
    test=:LRT,
    pvalfile = out_prefix * "_pval.csv",
    nullfile = out_prefix * ".log",
    geneticrowinds = sample_indices
)