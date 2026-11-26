# Exposure associations in single-cell data with exposomeSC

Abstract

exposomeSC relates exposures measured on donors to gene expression
within the annotated cell types of single-cell RNA-seq data. This
vignette goes through the package on two small real data sets: airway
epithelial cells of smokers and never smokers from the Human Lung Cell
Atlas, and blood cells of patients whose plasma was measured for
compounds of exogenous origin. It covers the container, the
donor-pseudobulk association model and its companions, and the analyses
built around it.

## Introduction

*[exposomeSC](https://bioconductor.org/packages/3.25/exposomeSC)*
relates exposures measured on donors to gene expression within the
annotated cell types of single-cell RNA-seq data. Its container, the
`SingleCellExposomeExperiment`, keeps a donor-level exposure table next
to the cells of a `SingleCellExperiment` and checks that the two agree.
Its main analysis, [`run_sc_exwas()`](../reference/run_sc_exwas.md),
sums the raw counts of each cell type within each donor and fits a model
in which the donors, not the cells, are the replicates.

Cells of one donor share that donor’s exposure, so they are not
independent observations of it; a test that treats them as such is too
confident ([Squair et al. 2021](#ref-squair2021)). The quantity that
[`run_sc_exwas()`](../reference/run_sc_exwas.md) estimates is the log2
fold change of a gene’s expression rate within a cell type per unit of
the exposure. A change in the number of cells of a type is a different
quantity and has a function of its own (Section @ref(composition)).

The package works with the classes of Bioconductor. The container
extends `SingleCellExperiment`, the pseudobulk counts are
`SummarizedExperiment` objects, and the results are `DataFrame`s; those
of [`run_sc_exwas()`](../reference/run_sc_exwas.md) carry the settings
of the run in their
[`metadata()`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html).
The models are fitted by
*[edgeR](https://bioconductor.org/packages/3.25/edgeR)*,
*[DESeq2](https://bioconductor.org/packages/3.25/DESeq2)* and
*[limma](https://bioconductor.org/packages/3.25/limma)*. In a workflow
the package comes after quality control and the annotation of cell
types, which it takes as given: its input is a `SingleCellExperiment`
with raw counts, a column that names the donor of each cell and a column
with the cell type.

Every association the package estimates is observational. Exposures such
as smoking are not assigned at random, and nothing in the package
identifies a causal effect.

The vignette uses the two data sets that the package ships. Both are
real and both are small on purpose: they show what the functions take
and what they return, and they are far too small to support conclusions.
The first, 352 airway epithelial cells of 12 donors, serves Sections
@ref(data) to @ref(composition). The second, 1,032 blood cells of 43
patients with measured plasma compounds and surface proteins, serves
Section @ref(incov).

Some sections use packages that exposomeSC only suggests:
*[lme4](https://CRAN.R-project.org/package=lme4)*,
*[lmerTest](https://CRAN.R-project.org/package=lmerTest)* and
*[pbkrtest](https://CRAN.R-project.org/package=pbkrtest)*,
*[mgcv](https://CRAN.R-project.org/package=mgcv)*,
*[qgcomp](https://CRAN.R-project.org/package=qgcomp)*,
*[metafor](https://CRAN.R-project.org/package=metafor)*,
*[ppcor](https://CRAN.R-project.org/package=ppcor)*,
*[variancePartition](https://bioconductor.org/packages/3.25/variancePartition)*
and *[fgsea](https://bioconductor.org/packages/3.25/fgsea)*. A chunk
that needs one of them is run when it is installed. Two further
vignettes cover the cross-omic networks and checks of a result by
permutation.

## Installation

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("exposomeSC")
```

## The airway data

``` r

library(exposomeSC)
library(SingleCellExperiment)
x <- readRDS(system.file("extdata", "hlca_airway_subset.rds",
    package = "exposomeSC"))
x
#> class: SingleCellExposomeExperiment 
#> dim: 20111 352 
#> metadata(1): source
#> assays(1): counts
#> rownames(20111): ENSG00000000003 ENSG00000000419 ... ENSG00000283103
#>   ENSG00000283125
#> rowData names(1): symbol
#> colnames(352): TCTTCCTTCAGAGCAG_T101 CTACTATCAACGGGTA_T101 ...
#>   CCTCTAGAGCAAGCCA_T90 CCTCATGGTCTACATG_T90
#> colData names(4): donor_id cell_type assay hlca_donor_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> exposureData(12 samples x 4 exposures): smoking, chemistry_v3, age, sex_male
#> sampleColumn: donor_id
```

The object holds the raw counts of 352 tracheal epithelial cells of the
Human Lung Cell Atlas core ([Sikkema et al. 2023](#ref-sikkema2023)),
from its study of smokers and never smokers ([Goldfarbmuren et al.
2020](#ref-goldfarbmuren2020)): six heavy smokers, with at least 15
pack-years, and six never smokers, with up to 10 basal, 10 secretory and
10 multiciliated-lineage cells per donor.
[`?hlca_airway_subset`](../reference/hlca_airway_subset.md) describes
the selection, and the script that made the file from the public H5AD is
installed with the package
(`system.file("scripts", package = "exposomeSC")`). The atlas is
distributed under the CC BY 4.0 licence, and the object carries its
provenance:

``` r

metadata(x)$source[c("data", "dataset_version", "licence", "references")]
#> $data
#> [1] "Human Lung Cell Atlas (HLCA) core, study Seibold_2020"
#> 
#> $dataset_version
#> [1] "688185ad-11c2-4172-a53a-f4f1f4076860"
#> 
#> $licence
#> [1] "CC BY 4.0"
#> 
#> $references
#> [1] "doi:10.1038/s41591-023-02327-2" "doi:10.1038/s41467-020-16239-z"
table(x$donor_id, x$cell_type)
#>       
#>        Basal Multiciliated lineage Secretory
#>   T101    10                    10        10
#>   T120    10                    10        10
#>   T126    10                    10        10
#>   T137    10                    10        10
#>   T153    10                    10        10
#>   T154    10                     2        10
#>   T164    10                    10        10
#>   T165    10                    10        10
#>   T166    10                    10        10
#>   T167    10                    10        10
#>   T85     10                    10        10
#>   T90     10                    10        10
```

The donor-level variables are the exposure data of the object, one row
per donor, with a description of each variable.
[`sampleColumn()`](../reference/SCEE-accessors.md) names the column of
[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
that assigns the cells to the donors.

``` r

exposureData(x)
#>      smoking chemistry_v3 age sex_male
#> T101       1            1  55        1
#> T120       1            0  57        0
#> T126       0            0  35        0
#> T137       0            0  27        0
#> T153       0            1  38        1
#> T154       1            1  61        1
#> T164       0            1  66        1
#> T165       0            1  64        1
#> T166       0            1  68        0
#> T167       1            1  66        0
#> T85        1            0  59        1
#> T90        1            1  44        1
exposureInfo(x)
#> DataFrame with 4 rows and 3 columns
#>       exposure        family                   unit
#>    <character>   <character>            <character>
#> 1      smoking tobacco smoke 1 = heavy smoker, 0 ..
#> 2 chemistry_v3     technical 1 = 10x 3' v3, 0 = 1..
#> 3          age   demographic                  years
#> 4     sex_male   demographic   1 = male, 0 = female
exposureVariables(x)
#> [1] "smoking"      "chemistry_v3" "age"          "sex_male"
sampleColumn(x)
#> [1] "donor_id"
head(cellSampleMap(x), 3)
#> DataFrame with 3 rows and 2 columns
#>                 cell_id   sample_id
#>             <character> <character>
#> 1 TCTTCCTTCAGAGCAG_T101        T101
#> 2 CTACTATCAACGGGTA_T101        T101
#> 3 ACAAGCTGTGCTGCAC_T101        T101
```

[`build_scee()`](../reference/build_scee.md) makes such an object from a
`SingleCellExperiment` and a numeric matrix with one row per donor. It
stops when a donor of the cells has no row in the matrix or a row has no
cells. Here the object is taken apart and built again:

``` r

sce <- SingleCellExperiment(
    assays = list(counts = counts(x)),
    colData = colData(x),
    rowData = rowData(x)
)
rebuilt <- build_scee(sce, exposureData(x), sample_col = "donor_id",
    exposure_info = exposureInfo(x))
identical(exposureData(rebuilt), exposureData(x))
#> [1] TRUE
```

Donor variables often arrive as columns of the cell annotations, with
the value of a donor repeated for each of its cells, as in the metadata
of a Seurat object.
[`seurat_to_exposure()`](../reference/seurat_to_exposure.md) reduces
such columns to the matrix that
[`build_scee()`](../reference/build_scee.md) takes:

``` r

cell_table <- data.frame(donor_id = x$donor_id,
    exposureData(x)[x$donor_id, ], row.names = colnames(x))
head(cell_table, 3)
#>                       donor_id smoking chemistry_v3 age sex_male
#> TCTTCCTTCAGAGCAG_T101     T101       1            1  55        1
#> CTACTATCAACGGGTA_T101     T101       1            1  55        1
#> ACAAGCTGTGCTGCAC_T101     T101       1            1  55        1
donor_table <- seurat_to_exposure(cell_table, sample_col = "donor_id",
    exposure_cols = exposureVariables(x))
identical(donor_table, exposureData(x))
#> [1] TRUE
```

[`donor_pseudobulk()`](../reference/donor_pseudobulk.md) returns the
counts that [`run_sc_exwas()`](../reference/run_sc_exwas.md) fits its
models to: one `SummarizedExperiment` per cell type, with one column per
donor and the exposure data as column data. Looking at them is the first
check of any result. The counts of *CYP1B1* in basal cells, a gene that
returns below, are:

``` r

pseudobulk <- donor_pseudobulk(x)
pseudobulk
#> List of length 3
#> names(3): Basal Multiciliated lineage Secretory
basal <- pseudobulk[["Basal"]]
cyp1b1 <- rownames(x)[rowData(x)$symbol == "CYP1B1"]
data.frame(
    smoking = basal$smoking,
    CYP1B1 = assay(basal, "counts")[cyp1b1, ],
    all_genes = colSums(assay(basal, "counts"))
)
#>      smoking CYP1B1 all_genes
#> T101       1      3    361506
#> T120       1     10    205282
#> T126       0      0    148021
#> T137       0      0    189713
#> T153       0      0    250864
#> T154       1     18    465165
#> T164       0      0    352726
#> T165       0      0    422411
#> T166       0      0    362814
#> T167       1     12    450820
#> T85        1      1    130760
#> T90        1      5    299829
```

## The association model

[`run_sc_exwas()`](../reference/run_sc_exwas.md) fits, for each cell
type, a donor-level model of the pseudobulk counts on the exposure and
the covariates. The default backend is the quasi-likelihood model of
*[edgeR](https://bioconductor.org/packages/3.25/edgeR)* ([Chen et al.
2025](#ref-chen2025)). Here the exposure is smoking, adjusted for the
10x chemistry, which differs between donors:

``` r

fit <- run_sc_exwas(x, exposure = "smoking", covariates = "chemistry_v3")
dim(fit)
#> [1] 32890    17
top <- head(fit, 8)
top$symbol <- rowData(x)[top$gene, "symbol"]
top[, c("symbol", "celltype", "log2FC", "statistic", "pvalue", "padj",
    "padj_global", "n_donors")]
#> DataFrame with 8 rows and 8 columns
#>        symbol              celltype    log2FC statistic      pvalue       padj
#>   <character>           <character> <numeric> <numeric>   <numeric>  <numeric>
#> 1       IL1R2                 Basal   2.96110   6.37915 5.17735e-07 0.00554391
#> 2        LDLR Multiciliated lineage   1.33962   6.01446 1.52521e-06 0.01105868
#> 3        KILH Multiciliated lineage   2.30558   5.79025 1.98558e-06 0.01105868
#> 4      CYP1B1                 Basal   5.97497   6.43537 2.86266e-06 0.01532670
#> 5      MIPOL1 Multiciliated lineage   2.05317   5.50351 5.04406e-06 0.01872859
#> 6        SAA4             Secretory  -3.46544  -5.41085 8.61588e-06 0.09514516
#> 7        PSCA Multiciliated lineage   1.32515   5.32518 1.02146e-05 0.02844509
#> 8      FGFBP1 Multiciliated lineage   2.81946   5.75107 1.33446e-05 0.02861956
#>   padj_global  n_donors
#>     <numeric> <integer>
#> 1   0.0170283        12
#> 2   0.0217686        11
#> 3   0.0217686        11
#> 4   0.0235382        12
#> 5   0.0331798        11
#> 6   0.0472294        12
#> 7   0.0479940        11
#> 8   0.0548628        11
```

The result has one row per tested gene and cell type, ordered by
`padj_global`. `padj` is the Benjamini-Hochberg adjustment within a cell
type and `padj_global` the adjustment across all tested
gene-by-cell-type pairs. For the edgeR backend `se` is `NA`: the
quasi-likelihood F test gives no standard error of the coefficient, and
`statistic` is the square root of the F statistic with the sign of the
log fold change. The settings and what was done in each cell type are
kept with the result:

``` r

parameters <- metadata(fit)$parameters
parameters[c("method", "se_method", "statistic_type")]
#> $method
#> [1] "edgeR"
#> 
#> $se_method
#> [1] "not_available_edgeR_QL"
#> 
#> $statistic_type
#> [1] "signed_sqrt_qlf"
diagnostics <- metadata(fit)$diagnostics
vapply(diagnostics, function(d) d$n_donors, numeric(1))
#>                 Basal Multiciliated lineage             Secretory 
#>                    12                    11                    12
vapply(diagnostics, function(d) d$n_genes_tested, numeric(1))
#>                 Basal Multiciliated lineage             Secretory 
#>                 10708                 11139                 11043
```

One donor has only two multiciliated-lineage cells, fewer than
`min_cells`, so that cell type is fitted with 11 donors.

These numbers must not be read as findings. With six donors per group
and at most ten cells per donor and cell type, a few counts decide a
result: the table of Section @ref(data) shows that *CYP1B1* has no count
in the basal cells of the never smokers and between 1 and 18 in those of
a smoker, and the size of such an estimate depends on little. The same
model on larger sets of cells of the same donors, which were analysed
when the data set was made, gave other lists of genes below a threshold.

## Other models for the same question

### Other backends

`method = "voom_dream"` fits the voom model of
*[limma](https://bioconductor.org/packages/3.25/limma)* through
*[variancePartition](https://bioconductor.org/packages/3.25/variancePartition)*
([Hoffman and Roussos 2021](#ref-hoffman2021)), and `method = "DESeq2"`
a Wald test of
*[DESeq2](https://bioconductor.org/packages/3.25/DESeq2)*. Both return a
standard error of the coefficient. The columns are the same for every
backend.

``` r

fit_voom <- run_sc_exwas(x, exposure = "smoking",
    covariates = "chemistry_v3", method = "voom_dream", celltypes = "Basal")
fit_voom$symbol <- rowData(x)[fit_voom$gene, "symbol"]
head(fit_voom[, c("symbol", "log2FC", "se", "statistic", "pvalue")], 3)
#> DataFrame with 3 rows and 5 columns
#>        symbol    log2FC        se statistic      pvalue
#>   <character> <numeric> <numeric> <numeric>   <numeric>
#> 1      CYP1B1   3.76235  0.477184   7.88450 3.54114e-08
#> 2       IL1R2   2.63220  0.472469   5.57117 9.20432e-06
#> 3         NTS  -1.40107  0.276195  -5.07276 3.26980e-05
```

### A negative binomial model with an offset

[`run_sc_exwas_pb_offset()`](../reference/run_sc_exwas_pb_offset.md)
fits a negative binomial model of each pseudobulk count with the
logarithm of the donor’s total count as offset ([Lee and Han
2024](#ref-lee2024)), for the genes named in `target_genes`, here three.
Its p-values are Wald tests and can be too small with few donors. It has
no `padj` within a cell type; `padj_global` is adjusted across the tests
that were made. It also has no `min_cells`: every donor with a cell of
the type enters the model, the donor with two multiciliated-lineage
cells included.

``` r

genes <- rownames(x)[match(c("ALDH3A1", "NQO1", "GPX2"), rowData(x)$symbol)]
offset_fit <- run_sc_exwas_pb_offset(x, exposure = "smoking",
    covariates = "chemistry_v3", target_genes = genes)
offset_fit$symbol <- rowData(x)[offset_fit$gene, "symbol"]
offset_fit[, c("symbol", "celltype", "log2FC", "se", "statistic", "pvalue",
    "padj_global")]
#> DataFrame with 9 rows and 7 columns
#>        symbol              celltype     log2FC        se statistic      pvalue
#>   <character>           <character>  <numeric> <numeric> <numeric>   <numeric>
#> 1     ALDH3A1 Multiciliated lineage  1.6639367  0.477171  3.487085 0.000488316
#> 2        GPX2 Multiciliated lineage  0.9285628  0.288305  3.220771 0.001278462
#> 3        NQO1 Multiciliated lineage  0.5981829  0.255481  2.341396 0.019211799
#> 4        GPX2             Secretory  0.6544349  0.415330  1.575699 0.115095177
#> 5        GPX2                 Basal -0.3673234  0.272701 -1.346983 0.177985593
#> 6     ALDH3A1             Secretory  0.5849155  0.676691  0.864376 0.387381420
#> 7        NQO1                 Basal -0.1472037  0.349979 -0.420607 0.674042118
#> 8        NQO1             Secretory  0.0879069  0.474085  0.185424 0.852896333
#> 9     ALDH3A1                 Basal  0.0800406  0.459262  0.174281 0.861644575
#>   padj_global
#>     <numeric>
#> 1  0.00439485
#> 2  0.00575308
#> 3  0.05763540
#> 4  0.25896415
#> 5  0.32037407
#> 6  0.58107213
#> 7  0.86164458
#> 8  0.86164458
#> 9  0.86164458
```

### A mixed model of the cells

[`run_sc_exwas_glmm()`](../reference/run_sc_exwas_glmm.md) keeps the
cells as observations. For each gene it fits a linear mixed model of the
log-normalised expression of the cells, with random intercepts for the
donor and for the donor within a cell type, and tests the effect of the
exposure in each cell type with denominator degrees of freedom by the
method of Satterthwaite, the default, or of Kenward and Roger, which its
help page advises for fewer than about 30 donors and which is used here.
It is much slower than the pseudobulk models, so it is meant for a few
genes.

``` r

mixed_fit <- run_sc_exwas_glmm(x, exposure = "smoking",
    covariates = "chemistry_v3", target_genes = genes,
    ddf = "Kenward-Roger")
mixed_fit$symbol <- rowData(x)[mixed_fit$gene, "symbol"]
mixed_fit[, c("symbol", "celltype", "log2FC", "se", "df", "pvalue",
    "padj_global")]
#> DataFrame with 9 rows and 7 columns
#>        symbol              celltype     log2FC        se        df    pvalue
#>   <character>           <character>  <numeric> <numeric> <numeric> <numeric>
#> 1     ALDH3A1 Multiciliated lineage  1.1937568  0.522799   13.4336 0.0392683
#> 2        NQO1 Multiciliated lineage  0.3922636  0.229207   23.1645 0.1003672
#> 3        GPX2 Multiciliated lineage  0.4269340  0.262276   17.3609 0.1215772
#> 4        GPX2                 Basal -0.3678811  0.256765   16.3480 0.1707738
#> 5        GPX2             Secretory  0.2085547  0.256765   16.3480 0.4283237
#> 6     ALDH3A1             Secretory  0.3727827  0.516562   12.9032 0.4833755
#> 7        NQO1                 Basal -0.1391103  0.223053   21.8187 0.5393143
#> 8     ALDH3A1                 Basal -0.1347923  0.516562   12.9032 0.7982561
#> 9        NQO1             Secretory  0.0514866  0.223053   21.8187 0.8196001
#>   padj_global
#>     <numeric>
#> 1    0.353415
#> 2    0.364732
#> 3    0.364732
#> 4    0.384241
#> 5    0.693404
#> 6    0.693404
#> 7    0.693404
#> 8    0.819600
#> 9    0.819600
```

The response is log2(counts per 10,000 + 1) of each cell, so the
coefficient, although the column is called `log2FC`, is not on the scale
of the pseudobulk estimates. `family = "nbinom"` fits a negative
binomial mixed model with a likelihood-ratio test for each gene and cell
type instead; that takes minutes for a few genes even on these data and
is not run here.

## More than one exposure, and differences between cell types

[`run_multi_exwas()`](../reference/run_multi_exwas.md) runs the model
for several exposures and adds `padj_all`, the adjustment across all
exposures and cell types. Age is used here as a second donor variable;
`celltypes` restricts the run to basal cells.

``` r

both <- run_multi_exwas(x, exposures = c("smoking", "age"),
    covariates = "chemistry_v3", celltypes = "Basal")
#> Running sc-ExWAS for: smoking
#> Running sc-ExWAS for: age
table(both$exposure)
#> 
#>     age smoking 
#>   10708   10708
both$symbol <- rowData(x)[both$gene, "symbol"]
head(both[, c("symbol", "exposure", "log2FC", "pvalue", "padj_global",
    "padj_all")], 3)
#> DataFrame with 3 rows and 6 columns
#>        symbol    exposure    log2FC      pvalue padj_global  padj_all
#>   <character> <character> <numeric>   <numeric>   <numeric> <numeric>
#> 1       IL1R2     smoking   2.96110 5.17735e-07  0.00554391 0.0110878
#> 2      CYP1B1     smoking   5.97497 2.86266e-06  0.01532670 0.0306534
#> 3      TSPAN8     smoking  -2.49851 2.03721e-05  0.05522154 0.1104431
```

A small p-value in one cell type and a large one in another do not show
that the effect differs between the two ([Gelman and Stern
2006](#ref-gelman2006)).
[`run_interaction_test()`](../reference/run_interaction_test.md) tests
the difference itself, for the genes in `target_genes` or, by default,
for the 50 most variable genes:

``` r

interaction <- run_interaction_test(x, exposure = "smoking",
    covariates = "chemistry_v3")
interaction$symbol <- rowData(x)[interaction$gene, "symbol"]
head(interaction[order(interaction$p_interaction),
    c("symbol", "p_interaction", "padj_interaction", "n_celltypes",
        "n_donors")], 3)
#> DataFrame with 3 rows and 5 columns
#>        symbol p_interaction padj_interaction n_celltypes  n_donors
#>   <character>     <numeric>        <numeric>   <integer> <integer>
#> 1      BPIFA1     0.0142263         0.711314           3        12
#> 2       WDR38     0.0675888         0.957566           3        12
#> 3     DYNLRB2     0.0972912         0.957566           3        12
```

## Gene sets

[`run_gsea()`](../reference/run_gsea.md) ranks the genes of each cell
type by the test statistic and runs
*[fgsea](https://bioconductor.org/packages/3.25/fgsea)* on the ranking.
Gene sets are named lists of gene identifiers as in the row names of the
object. The three sets below are short lists written for this example,
not a curated collection; they are mapped from symbols to the Ensembl
identifiers of the data.

``` r

symbols <- list(
    xenobiotic_response = c("CYP1B1", "CYP1A1", "AHRR", "ALDH3A1", "NQO1",
        "AKR1C1", "AKR1C2", "AKR1B10", "GPX2", "TXNRD1", "CBR1", "UGT1A6"),
    cilia = c("FOXJ1", "DNAH5", "DNAI1", "CCDC39", "CCDC40", "RSPH1",
        "TEKT1", "DNAH9", "CFAP43", "SPAG6"),
    secretory = c("MUC5AC", "MUC5B", "SCGB1A1", "TFF3", "BPIFB1", "AGR2",
        "SPDEF", "FOXA3")
)
gene_sets <- lapply(symbols, function(set) {
    rownames(x)[rowData(x)$symbol %in% set]
})
lengths(gene_sets)
#> xenobiotic_response               cilia           secretory 
#>                  12                   9                   8
```

``` r

set.seed(1)
enrichment <- run_gsea(fit, gene_sets = gene_sets, min_size = 5L)
enrichment[, c("pathway", "celltype", "NES", "pval", "padj_global", "size")]
#> DataFrame with 7 rows and 6 columns
#>               pathway              celltype       NES        pval padj_global
#>           <character>           <character> <numeric>   <numeric>   <numeric>
#> 1           secretory                 Basal -0.827794 0.697318008  0.69731801
#> 2 xenobiotic_response                 Basal  1.158095 0.284536082  0.33195876
#> 3               cilia Multiciliated lineage -1.564555 0.045194942  0.07909115
#> 4           secretory Multiciliated lineage  1.838749 0.002449421  0.00857297
#> 5 xenobiotic_response Multiciliated lineage  1.927211 0.000589535  0.00412674
#> 6           secretory             Secretory -1.640713 0.021349762  0.04981611
#> 7 xenobiotic_response             Secretory  1.278132 0.173410405  0.24277457
#>        size
#>   <integer>
#> 1         7
#> 2        10
#> 3         9
#> 4         7
#> 5         8
#> 6         8
#> 7         9
```

`size` is the number of genes of a set that were tested in the cell
type. A set with fewer than `min_size` of them is left out for that cell
type.

## Combining cohorts

[`run_meta_exwas()`](../reference/run_meta_exwas.md) combines the
estimates of independent cohorts by inverse-variance weighting. It needs
results with standard errors, so results of the edgeR backend are
refused, and its random-effects synthesis (`method = "random"`) needs at
least three cohorts. The airway data are one study; to show the
function, the four donors measured with the 10x v2 chemistry and the
eight measured with v3 are treated as two cohorts, each with too few
donors for an analysis of its own, and the genes of the sets above are
combined.

``` r

v3 <- rownames(exposureData(x))[exposureData(x)[, "chemistry_v3"] == 1]
cohorts <- list(
    v2 = x[, !x$donor_id %in% v3],
    v3 = x[, x$donor_id %in% v3]
)
cohort_fits <- lapply(cohorts, function(cohort) {
    cohort_fit <- run_sc_exwas(cohort, exposure = "smoking",
        method = "voom_dream", celltypes = "Basal", min_donors = 4L)
    cohort_fit[cohort_fit$gene %in% unlist(gene_sets), ]
})
meta <- run_meta_exwas(cohort_fits, cohort_names = names(cohort_fits))
meta$symbol <- rowData(x)[meta$gene, "symbol"]
head(meta[order(meta$meta_pvalue), c("symbol", "meta_effect", "meta_se",
    "meta_pvalue", "n_cohorts", "direction")], 3)
#> DataFrame with 3 rows and 6 columns
#>        symbol meta_effect   meta_se meta_pvalue n_cohorts   direction
#>   <character>   <numeric> <numeric>   <numeric> <integer> <character>
#> 1      CYP1B1     3.76661  0.493605 2.33272e-14         2          ++
#> 2     SCGB1A1    -3.01158  0.532967 1.59893e-08         2          --
#> 3      AKR1C1    -1.01850  0.315096 1.22781e-03         2          --
```

Subsetting the columns of the object keeps the exposure data of the
remaining donors. The combination treats the standard error of each
cohort as known. With two donors against two in one cohort and four
against four in the other they are anything but, so the p-values of this
table show how the function is called and nothing about the genes.

## Cell composition

[`run_exposure_composition()`](../reference/run_exposure_composition.md)
asks whether the proportions of the cell types differ with the exposure,
with the donors as the units. The cells of the airway data were drawn in
equal numbers per donor and cell type, so they carry no information
about composition. The provenance record of the object holds the numbers
of cells of each donor and cell type in the atlas from which they were
drawn, and an object with one label per cell and no counts is enough for
the function:

``` r

eligible <- metadata(x)$source$eligible_cells
eligible
#>      Basal Multiciliated lineage Secretory
#> T101  2957                   140       112
#> T120  1342                   185       216
#> T126  1454                   192       100
#> T137   911                   282       119
#> T153  2690                   251       394
#> T154  1413                     2       210
#> T164  2913                    79       194
#> T165  1983                    51        33
#> T166  2079                   299       191
#> T167  2220                   163       386
#> T85   3053                   108       316
#> T90   1873                    90       216
cells <- paste0("cell", seq_len(sum(eligible)))
labels <- DataFrame(
    donor_id = rep(rownames(eligible)[row(eligible)], eligible),
    cell_type = rep(colnames(eligible)[col(eligible)], eligible),
    row.names = cells
)
no_counts <- matrix(0, 0, length(cells), dimnames = list(NULL, cells))
atlas_cells <- build_scee(
    SingleCellExperiment(assays = list(counts = no_counts), colData = labels),
    exposureData(x), sample_col = "donor_id"
)
composition <- run_exposure_composition(atlas_cells, exposure = "smoking",
    covariates = "chemistry_v3")
composition[, c("log2_ratio_change", "se", "pvalue", "padj",
    "median_observed_proportion")]
#> DataFrame with 3 rows and 5 columns
#>                       log2_ratio_change        se    pvalue      padj
#>                               <numeric> <numeric> <numeric> <numeric>
#> Basal                          0.379082  0.710338  0.606501  0.606501
#> Multiciliated lineage         -1.710891  1.128324  0.163747  0.359128
#> Secretory                      1.331809  0.774701  0.119709  0.359128
#>                       median_observed_proportion
#>                                        <numeric>
#> Basal                                  0.8461646
#> Multiciliated lineage                  0.0512467
#> Secretory                              0.0907921
```

`log2_ratio_change` is the coefficient of the exposure for the log2
ratio of a cell type to the geometric mean of the others, with a
heteroscedasticity-consistent standard error, and `padj` is the Holm
adjustment over the cell types. A test of the composition as a whole is
in `metadata(composition)$omnibus`. The ratios are among these three
cell types as they were captured; they say nothing about absolute
numbers of cells in the tissue.

## A cohort with measured exposures and proteins

The airway data have a smoking status and no measurement of what the
donors were exposed to. The second data set of the package has such
measurements: 1,032 blood cells of 43 patients of the INCOV cohort ([Su
et al. 2020](#ref-su2020)), 14 monocytes and 10 T cells each, with the
counts of 192 antibodies for the same cells, four of them isotype
controls, and the plasma levels of 13 compounds of exogenous origin
among the donor variables. All donors had COVID-19 and were not in
hospital when the blood was drawn.
[`?incov_pbmc_subset`](../reference/incov_pbmc_subset.md) gives the
sources, their licences and every step from the public files to the
object.

``` r

y <- readRDS(system.file("extdata", "incov_pbmc_subset.rds",
    package = "exposomeSC"))
y
#> class: SingleCellExposomeExperiment 
#> dim: 15406 1032 
#> metadata(3): source plasma_metabolites metabolite_info
#> assays(1): counts
#> rownames(15406): A1BG A1BG-AS1 ... ZYX ZZEF1
#> rowData names(1): symbol
#> colnames(1032): INCOV039_AAGACCTAGGGTCTCC INCOV039_ACTGTCCGTGAGGCTA ...
#>   INCOV135_GGCGTGTGTACCGTAT INCOV135_GTATCTTAGAAACGCC
#> colData names(4): donor_id cell_type sample_id archive_cell
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(1): CITE
#> exposureData(43 samples x 16 exposures): ever_smoker, age, sex_male, cotinine, caffeine ...
#> sampleColumn: donor_id
table(y$cell_type)
#> 
#> Monocyte   T cell 
#>      602      430
as.data.frame(exposureInfo(y)[, c("exposure", "family", "source_name")])
#>                 exposure         family                      source_name
#> 1            ever_smoker  tobacco smoke                Cigarette Smoking
#> 2                    age    demographic                              Age
#> 3               sex_male    demographic                              Sex
#> 4               cotinine        tobacco                         cotinine
#> 5               caffeine       caffeine                         caffeine
#> 6           paraxanthine       caffeine                     paraxanthine
#> 7            theobromine       caffeine                      theobromine
#> 8                   pfos perfluoroalkyl  perfluorooctanesulfonate (PFOS)
#> 9                   pfoa perfluoroalkyl        perfluorooctanoate (PFOA)
#> 10 methylparaben_sulfate        paraben methyl-4-hydroxybenzoate sulfate
#> 11 propylparaben_sulfate        paraben propyl 4-hydroxybenzoate sulfate
#> 12 hydroxychlorothalonil      pesticide          4-hydroxychlorothalonil
#> 13      naphthol_sulfate    naphthalene               2-naphthol sulfate
#> 14     ethyl_glucuronide        alcohol                ethyl glucuronide
#> 15             saccharin      sweetener                        saccharin
#> 16            acesulfame      sweetener                       acesulfame
```

### Compounds without a value

The plasma levels are relative abundances from untargeted mass
spectrometry, scaled by the authors so that the median of a compound
within a batch is 1. For some donors the table has no value of a
compound:

``` r

colSums(is.na(exposureData(y)[, c("caffeine", "paraxanthine", "theobromine",
    "pfos", "pfoa", "saccharin")]))
#>     caffeine paraxanthine  theobromine         pfos         pfoa    saccharin 
#>            3            3            0            1            1           16
```

The table does not say why a value is missing. In such measurements the
usual reason is that the compound was not detected, but it is not the
only one, and the data can speak against it. Paraxanthine is the main
metabolite of caffeine:

``` r

pair <- exposureData(y)[, c("caffeine", "paraxanthine")]
round(pair[!complete.cases(pair), ], 2)
#>          caffeine paraxanthine
#> INCOV039     9.82           NA
#> INCOV068     4.24           NA
#> INCOV081     2.87           NA
#> INCOV109       NA         0.04
#> INCOV116       NA         0.05
#> INCOV134       NA         0.27
```

The donors without a value of caffeine have little paraxanthine, which
fits caffeine that was not detected. The donors without a value of
paraxanthine have caffeine well above the median of their batch, which
is 1; their paraxanthine is unlikely to have been below detection, so
paraxanthine is left out of what follows.

[`exposure_impute_lod()`](../reference/exposure_impute_lod.md) replaces
values below a limit of detection (LOD), by substitution, by regression
on order statistics (`"ros"`) or by multiple imputation (`"multiple"`).
`censored` marks the values that were below the limit; a missing value
that is not marked stays missing.

Three things are assumed here to show the function, and none is in the
source. The section takes every missing value of the five compounds
below to be a compound that was not detected. The table reports no
limits of detection; the smallest value of a compound among all samples
of the study, kept in `exposureInfo(y)$lowest_reported`, stands in for
the limit, although the values were scaled batch by batch and one limit
for all batches is itself an approximation. And `"ros"` takes the
logarithms of the levels to be normally distributed.

``` r

compounds <- c("caffeine", "theobromine", "pfos", "pfoa", "saccharin")
reported <- exposureData(y)[, compounds]
```

``` r

info <- exposureInfo(y)
lowest <- setNames(info$lowest_reported, info$exposure)[compounds]
signif(lowest, 2)
#>    caffeine theobromine        pfos        pfoa   saccharin 
#>      0.0190      0.0016      0.0720      0.0650      0.0910
filled <- exposure_impute_lod(reported, lod = lowest, method = "ros",
    censored = is.na(reported))
#> caffeine: 3 below LOD, 0 missing, 40 detected
#> theobromine: 0 below LOD, 0 missing, 43 detected
#> pfos: 1 below LOD, 0 missing, 42 detected
#> pfoa: 1 below LOD, 0 missing, 42 detected
#> saccharin: 16 below LOD, 0 missing, 27 detected
attr(filled, "lod_summary")
#>      exposure n_detected n_censored n_missing method
#> 1    caffeine         40          3         0    ros
#> 2 theobromine         43          0         0   none
#> 3        pfos         42          1         0    ros
#> 4        pfoa         42          1         0    ros
#> 5   saccharin         27         16         0    ros
fill_over_limit <- function(compound) {
    gone <- is.na(reported[, compound])
    round(sort(filled[gone, compound]) / lowest[[compound]], 3)
}
fill_over_limit("caffeine")
#> INCOV109 INCOV116 INCOV134 
#>    0.999    0.999    0.999
fill_over_limit("saccharin")
#> INCOV039 INCOV042 INCOV057 INCOV058 INCOV078 INCOV081 INCOV088 INCOV094 
#>    0.047    0.100    0.153    0.210    0.271    0.337    0.408    0.485 
#> INCOV100 INCOV102 INCOV103 INCOV107 INCOV116 INCOV131 INCOV132 INCOV134 
#>    0.569    0.661    0.760    0.869    0.987    0.999    0.999    0.999
```

`"ros"` gives one filled-in value for each missing one. The fills follow
the distribution fitted to the reported values and say nothing about
which donor had the larger level, which suits a summary and not an
association. With a limit as low as the one assumed here, the three
fills of caffeine all lie just below it, and only saccharin, with 16
missing values, gets a spread of them. For an association analysis
[`?exposure_impute_lod`](../reference/exposure_impute_lod.md) advises
`method = "multiple"`, which returns several filled-in data sets to
analyse one by one and combine; the sections below use the single
`"ros"` fill all the same, to stay short.

The logarithms of the filled levels replace the reported levels in the
object; `exposureData<-` checks that the donors still agree with the
cells, and `exposureInfo<-` records what the values now are.

``` r

levels <- exposureData(y)
levels[, compounds] <- log2(filled)
exposureData(y) <- levels
round(head(exposureData(y)[, compounds], 3), 2)
#>          caffeine theobromine  pfos  pfoa saccharin
#> INCOV039     3.30       -0.19 -0.92 -0.78     -7.88
#> INCOV042    -0.03        0.67  1.35 -1.22     -6.79
#> INCOV056     2.66       -0.06 -1.06 -0.07      1.74
info$unit[match(compounds, info$exposure)] <-
    "log2 of the relative abundance, missing values filled in"
exposureInfo(y) <- info
unique(exposureInfo(y)$unit)
#> [1] "1 = current or former smoker, 0 = never smoker"          
#> [2] "years"                                                   
#> [3] "1 = male, 0 = female"                                    
#> [4] "relative abundance, median of a batch = 1"               
#> [5] "log2 of the relative abundance, missing values filled in"
```

### Dose and response

[`run_dose_response()`](../reference/run_dose_response.md) compares
polynomial models of a gene’s pseudobulk expression on a continuous
exposure within one cell type, and
[`run_dose_response_gam()`](../reference/run_dose_response_gam.md) fits
a smooth curve with *[mgcv](https://CRAN.R-project.org/package=mgcv)*.
The dose here is the plasma level of caffeine, and the genes are 12
interferon-stimulated genes with counts in most donors, which serve as
the outcome in the rest of this section because the donors had a viral
infection.

``` r

interferon <- c("IFI6", "ISG15", "MX1", "IFITM3", "STAT1", "XAF1", "EPSTI1",
    "LY6E", "IRF7", "ISG20", "MX2", "EIF2AK2")
dose <- run_dose_response(y, exposure = "caffeine", celltype = "Monocyte",
    target_genes = interferon)
head(dose[order(dose$p_nonlinear), c("gene", "best_model", "linear_coef",
    "p_nonlinear", "n_donors")], 4)
#> DataFrame with 4 rows and 5 columns
#>          gene  best_model linear_coef p_nonlinear  n_donors
#>   <character> <character>   <numeric>   <numeric> <integer>
#> 1       ISG15   quadratic  -0.1056280    0.139063        43
#> 2      EPSTI1      linear  -0.1623649    0.369017        43
#> 3     EIF2AK2      linear  -0.0635761    0.373550        43
#> 4        IFI6      linear  -0.2664277    0.412985        43
```

`best_model` is the polynomial with the lowest AIC, and `p_nonlinear` an
F test of the cubic against the linear model; it is not adjusted for the
number of genes.

``` r

smooth <- run_dose_response_gam(y, exposure = "caffeine",
    celltype = "Monocyte", target_genes = interferon[1:4])
smooth[, c("gene", "edf", "p_smooth", "p_nonlinear", "R2_loocv_gam")]
#> DataFrame with 4 rows and 5 columns
#>          gene       edf  p_smooth p_nonlinear R2_loocv_gam
#>   <character> <numeric> <numeric>   <numeric>    <numeric>
#> 1        IFI6   1.00011 0.0385404   0.4335986   0.00947767
#> 2       ISG15   1.93161 0.1655689   0.0769601  -0.01951944
#> 3         MX1   1.00006 0.5223509   0.7782516  -0.07199001
#> 4      IFITM3   1.00002 0.0900116   0.4186522  -0.07171163
```

### Mixtures

[`run_sc_mixture()`](../reference/run_sc_mixture.md) scores each
exposure by its quantile and regresses the pseudobulk expression of each
target gene on the scores; it returns the sum of the coefficients for
each gene and, in the metadata, the relative size of the coefficients of
each exposure. It is a screen, without intervals or tests.
[`run_mixture_qgcomp()`](../reference/run_mixture_qgcomp.md) fits
quantile g-computation ([Keil et al. 2020](#ref-keil2020)) with
*[qgcomp](https://CRAN.R-project.org/package=qgcomp)* to the mean
expression of the target genes and returns an interval and a p-value for
the effect of raising all exposures by one quantile.

``` r

mixture_of <- c("caffeine", "theobromine", "pfos", "pfoa")
mixture <- run_sc_mixture(y, exposures = mixture_of, celltype = "Monocyte",
    target_genes = interferon, covariates = c("age", "sex_male"))
head(mixture[order(-abs(mixture$mixture_effect)), ], 3)
#> DataFrame with 3 rows and 2 columns
#>          gene mixture_effect
#>   <character>      <numeric>
#> 1        IFI6      -0.987597
#> 2      EPSTI1      -0.635306
#> 3      IFITM3      -0.528017
round(metadata(mixture)$weights, 2)
#>    caffeine theobromine        pfos        pfoa 
#>        0.30        0.23        0.16        0.31
```

``` r

joint <- run_mixture_qgcomp(y, exposures = mixture_of, celltype = "Monocyte",
    target_genes = interferon, covariates = c("age", "sex_male"))
joint[, c("mixture_effect", "ci_lower", "ci_upper", "mixture_pvalue",
    "n_donors")]
#> DataFrame with 1 row and 5 columns
#>   mixture_effect  ci_lower  ci_upper mixture_pvalue  n_donors
#>        <numeric> <numeric> <numeric>      <numeric> <integer>
#> 1      -0.328523  -1.02009   0.36304       0.341756        43
```

### Genes and proteins

[`run_cell_coupling()`](../reference/run_cell_coupling.md) is for data
with surface proteins measured in the same cells, kept as an alternative
experiment. Within each donor it computes the Spearman partial
correlation between a gene and a protein across the cells of one cell
type, given the RNA and protein library sizes; a meta-regression across
the donors then asks whether this correlation differs with the exposure.
It takes every gene with every protein; the rows of five genes with
their own proteins are shown.

``` r

pairs <- data.frame(
    gene = c("HLA-DRA", "CD36", "ITGAM", "FCGR3A", "CD86"),
    protein = c("HLA-DR", "CD36-1", "CD11B", "CD16", "CD86-1")
)
coupling <- run_cell_coupling(y, genes = pairs$gene,
    proteins = pairs$protein, exposure = "ever_smoker", celltype = "Monocyte",
    min_cells = 14L)
#> [exposomeSC] Cell coupling: Monocyte, 5 genes x 5 proteins = 25 pairs, 2 confounders
#> [exposomeSC] 25 pairs tested for Monocyte
own <- paste(coupling$gene, coupling$protein) %in%
    paste(pairs$gene, pairs$protein)
coupling[own, c("gene", "protein", "n_donors", "mean_r", "beta1", "se1",
    "pval1")]
#> DataFrame with 5 rows and 7 columns
#>          gene     protein  n_donors    mean_r      beta1       se1     pval1
#>   <character> <character> <integer> <numeric>  <numeric> <numeric> <numeric>
#> 1     HLA-DRA      HLA-DR        43 0.2696685 0.32435086  0.130771 0.0173242
#> 2        CD36      CD36-1        43 0.1099379 0.11703345  0.127730 0.3648920
#> 3       ITGAM       CD11B        39 0.1128381 0.21469605  0.122741 0.0885514
#> 4      FCGR3A        CD16        39 0.3816497 0.10903979  0.131641 0.4128062
#> 5        CD86      CD86-1        39 0.0581079 0.00191838  0.119037 0.9872286
```

`mean_r` is the mean of the partial correlation over the donors, and
`beta1` the slope, across donors, of its Fisher z on the exposure. The
default `min_cells` is 30 cells per donor; the data have 14 monocytes
per donor, too few for an analysis that is meant to be read.

## Related packages

exposomeSC adds a container for donor-level exposures and analyses that
are specific to exposure studies; for the model fits themselves it calls
other packages. The pseudobulk counts are summed by
*[scrapper](https://bioconductor.org/packages/3.25/scrapper)*, and the
models are those of
*[edgeR](https://bioconductor.org/packages/3.25/edgeR)*,
*[DESeq2](https://bioconductor.org/packages/3.25/DESeq2)* and
*[limma](https://bioconductor.org/packages/3.25/limma)* with
*[variancePartition](https://bioconductor.org/packages/3.25/variancePartition)*.
The packages below cover neighbouring ground.

- *[muscat](https://bioconductor.org/packages/3.25/muscat)* analyses
  differential states in multi-sample, multi-group single-cell data,
  with pseudobulk and mixed-model methods.
  [`run_sc_exwas()`](../reference/run_sc_exwas.md) fits a pseudobulk
  model of the same kind for a donor-level exposure that is kept in the
  object.
- *[dreamlet](https://bioconductor.org/packages/3.25/dreamlet)* tests
  pseudobulk differential expression in complex study designs, with the
  mixed models of variancePartition. `method = "voom_dream"` calls
  variancePartition itself, without a random effect; it does not use
  dreamlet.
- *[rexposome](https://bioconductor.org/packages/3.25/rexposome)* and
  *[omicRexposome](https://bioconductor.org/packages/3.25/omicRexposome)*
  hold exposome data and associate them with bulk omic data.
  [`as_scee()`](../reference/as_scee.md) takes the exposures of an
  `ExposomeSet` of rexposome into the container.
- *[speckle](https://bioconductor.org/packages/3.25/speckle)*,
  *[crumblr](https://bioconductor.org/packages/3.25/crumblr)* and
  *[sccomp](https://bioconductor.org/packages/3.25/sccomp)* test
  differences in the composition of cell types, as
  [`run_exposure_composition()`](../reference/run_exposure_composition.md)
  does; its help page describes how it differs from each.

A study whose question is a comparison of groups of samples, without
exposure measurements to keep with the cells, is served by the first two
packages directly.

## Session information

``` r

sessionInfo()
#> R Under development (unstable) (2026-11-25 r90991)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SingleCellExperiment_1.36.0 SummarizedExperiment_1.44.0
#>  [3] Biobase_2.74.0              GenomicRanges_1.66.0       
#>  [5] Seqinfo_1.4.0               IRanges_2.48.0             
#>  [7] S4Vectors_0.52.0            BiocGenerics_0.60.0        
#>  [9] generics_0.1.4              MatrixGenerics_1.26.0      
#> [11] matrixStats_1.5.0           exposomeSC_0.99.4          
#> [13] BiocStyle_2.42.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] mathjaxr_2.0-0           RColorBrewer_1.1-3       jsonlite_2.0.0          
#>   [4] magrittr_2.0.5           farver_2.1.2             AER_1.2-17              
#>   [7] nloptr_2.2.1             rmarkdown_2.32           fs_2.1.0                
#>  [10] ragg_1.5.2               vctrs_0.7.3              minqa_1.2.8             
#>  [13] htmltools_0.5.9          S4Arrays_1.14.0          BiocNeighbors_2.8.0     
#>  [16] broom_1.0.13             SparseArray_1.13.5       Formula_1.2-6           
#>  [19] ppcor_1.1                variancePartition_1.44.0 sass_0.5.0              
#>  [22] KernSmooth_2.23-27       bslib_0.12.0             desc_1.4.3              
#>  [25] pbkrtest_0.5.5           plyr_1.8.9               sandwich_3.1-3          
#>  [28] rootSolve_1.8.2.4        zoo_1.9-1                cachem_1.1.0            
#>  [31] lifecycle_1.0.5          iterators_1.0.14         pkgconfig_2.0.3         
#>  [34] Matrix_1.7-6             R6_2.6.1                 fastmap_1.2.0           
#>  [37] rbibutils_2.4.1          digest_0.6.39            numDeriv_2016.8-1.1     
#>  [40] textshaping_1.1.0        beachmat_2.30.0          metadat_1.6-0           
#>  [43] abind_1.4-8              mgcv_1.9-5               compiler_4.7.0          
#>  [46] aod_1.3.3                S7_1.0.0                 backports_1.5.1         
#>  [49] metafor_5.2-1            BiocParallel_1.48.0      carData_3.0-6           
#>  [52] gplots_3.3.0             MASS_7.3-66              DelayedArray_0.40.0     
#>  [55] corpcor_1.6.10           gtools_3.9.5             caTools_1.18.4          
#>  [58] tools_4.7.0              lmtest_0.9-40            otel_0.3.0              
#>  [61] nnet_7.3-21              remaCor_0.0.20           glue_1.8.1              
#>  [64] nlme_3.1-172             grid_4.7.0               reshape2_1.4.5          
#>  [67] fgsea_1.40.0             gtable_0.3.6             tidyr_1.3.2             
#>  [70] data.table_1.18.6.1      car_3.1-5                XVector_0.54.0          
#>  [73] pillar_1.11.1            stringr_1.6.0            limma_4.0.0             
#>  [76] splines_4.7.0            dplyr_1.2.1              lattice_0.23-1          
#>  [79] survival_3.8-12          tidyselect_1.2.1         qgcomp_2.19.6           
#>  [82] knitr_1.52               reformulas_0.4.4         gridExtra_2.3.1         
#>  [85] scrapper_1.8.0           bookdown_0.48            edgeR_5.0.0             
#>  [88] RhpcBLASctl_0.23-42      xfun_0.62                statmod_1.5.2           
#>  [91] arm_1.15-4               stringi_1.8.9            yaml_2.3.12             
#>  [94] boot_1.3-32              evaluate_1.0.5           codetools_0.2-20        
#>  [97] tibble_3.3.1             BiocManager_1.30.27      cli_3.6.6               
#> [100] systemfonts_1.3.4        Rdpack_2.6.6             jquerylib_0.1.4         
#> [103] pscl_1.5.9               Rcpp_1.1.2               EnvStats_3.1.0          
#> [106] coda_0.19-4.1            parallel_4.7.0           pkgdown_2.2.1           
#> [109] ggplot2_4.0.4            bitops_1.1-0             lme4_2.1-0              
#> [112] mvtnorm_1.5-0            lmerTest_3.2-1           scales_1.4.0            
#> [115] purrr_1.2.2              fANCOVA_0.6-1            rlang_1.3.0             
#> [118] cowplot_1.2.0            fastmatch_1.1-8
```

## References

Chen, Yunshun, Lizhong Chen, Aaron T. L. Lun, Pedro L. Baldoni, and
Gordon K. Smyth. 2025. “edgeR V4: Powerful Differential Analysis of
Sequencing Data with Expanded Functionality and Improved Support for
Small Counts and Larger Datasets.” *Nucleic Acids Research* 53 (2):
gkaf018. <https://doi.org/10.1093/nar/gkaf018>.

Gelman, Andrew, and Hal Stern. 2006. “The Difference Between
‘Significant’ and ‘Not Significant’ Is Not Itself Statistically
Significant.” *The American Statistician* 60 (4): 328–31.
<https://doi.org/10.1198/000313006X152649>.

Goldfarbmuren, Katherine C. et al. 2020. “Dissecting the Cellular
Specificity of Smoking Effects and Reconstructing Lineages in the Human
Airway Epithelium.” *Nature Communications* 11: 2485.
<https://doi.org/10.1038/s41467-020-16239-z>.

Hoffman, Gabriel E., and Panos Roussos. 2021. “Dream: Powerful
Differential Expression Analysis for Repeated Measures Designs.”
*Bioinformatics* 37 (2): 192–201.
<https://doi.org/10.1093/bioinformatics/btaa687>.

Keil, Alexander P., Jessie P. Buckley, Katie M. O’Brien, Kelly K.
Ferguson, Shanshan Zhao, and Alexandra J. White. 2020. “A Quantile-Based
g-Computation Approach to Addressing the Effects of Exposure Mixtures.”
*Environmental Health Perspectives* 128 (4): 047004.
<https://doi.org/10.1289/EHP5838>.

Lee, Hanbin, and Buhm Han. 2024. “Pseudobulk with Proper Offsets Has the
Same Statistical Properties as Generalized Linear Mixed Models in
Single-Cell Case-Control Studies.” *Bioinformatics* 40 (8): btae498.
<https://doi.org/10.1093/bioinformatics/btae498>.

Sikkema, Lisa et al. 2023. “An Integrated Cell Atlas of the Lung in
Health and Disease.” *Nature Medicine* 29 (6): 1563–77.
<https://doi.org/10.1038/s41591-023-02327-2>.

Squair, Jordan W. et al. 2021. “Confronting False Discoveries in
Single-Cell Differential Expression.” *Nature Communications* 12: 5692.
<https://doi.org/10.1038/s41467-021-25960-2>.

Su, Yapeng, Daniel Chen, Dan Yuan, et al. 2020. “Multi-Omics Resolves a
Sharp Disease-State Shift Between Mild and Moderate COVID-19.” *Cell*
183 (6): 1479–1495.e20. <https://doi.org/10.1016/j.cell.2020.10.037>.
