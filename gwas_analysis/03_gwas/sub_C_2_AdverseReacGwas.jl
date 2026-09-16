"""
# Sub routine script for OrdinalGWAS
"""

using OrdinalGWAS
using CSV, DataFrames
using StatsModels

function ordinal_lrt(;
    phenotype::AbstractString,
    cov_list::AbstractString,
    ord_file::AbstractString,
    bim_name::AbstractString,
    out_prefix::AbstractString,
    snp_list::Union{Nothing, Integer, AbstractString} = nothing
    )

    try
        # Prepare Environment variables
        fam_file = bim_name * ".fam"

        # Determine delimiter
        first_line = readline(fam_file)
        delim = occursin("\t", first_line) ? "\t" : " "

        cov_df = CSV.read(ord_file, DataFrame)
        fam_df = CSV.read(fam_file, DataFrame, header=false, delim=delim)
        plink_sample_ids = fam_df[:, 2]

        # Extract index in plink file from ord_file
        sample_ids = cov_df[:, "iid"]
        sample_indices = findall(in(sample_ids), plink_sample_ids)

        # Define formula
        response = Term(Symbol(phenotype))
        covariates = [Term(Symbol(cov)) for cov in split(cov_list, "+")]

        formula_expr = response ~ sum(covariates)
        # formula_expr = @eval @formula($(Symbol(phenotype)) ~ $(Meta.parse(cov_list)))

        ordinalgwas(
            formula_expr,
            ord_file,
            bim_name,
            test=:LRT,
            pvalfile = out_prefix * "_pval.csv",
            nullfile = out_prefix * ".log",
            geneticrowinds = sample_indices
        )

    catch e
        error_log = out_prefix * ".error"

        println("Error occurred: ", e)
        println("Stacktrace: ")
        Base.show_backtrace(stderr, catch_backtrace())

        open(error_log, "a") do io
            println(io, "Error occurred: ", e)
            println(io, "Stacktrace: ")
            Base.show_backtrace(io, catch_backtrace())
        end
end

end

precompile(
    ordinal_lrt, (String, String, String, String, String, Union{Nothing, Integer, String})
)