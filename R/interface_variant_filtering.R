createFilteringUI <- function() {
    tagList(
        h2("Filtering Parameters"),
        column(
            width = 6,
            style = paste0(
                "background-color: #f7f7f7; ",
                "padding: 10px; margin-bottom: 20px;"),
            createNumericInputWithPopover("depth_threshold", "Depth Threshold:",
                "Enter the minimum depth of coverage for filtering.",
                value = 10, min = 1),
            createNumericInputWithPopover("genotype_quality_threshold",
                "Genotype Quality Threshold:",
                "Enter the minimum genotype quality score.",
                value = 30, min = 1),
            createNumericInputWithPopover("vaf_ref", "VAF Reference Threshold:",
                "Define the upper VAF threshold for wildtype variants",
                value = 5, min = 0, max = 100, step = 10),
            createNumericInputWithPopover("vaf_hom",
                "VAF Homozygous Threshold:",
                "Define the lower VAF threshold for homozygous variants.",
                value = 95, min = 0, max = 100, step = 10)),
        column(
            width = 6,
            style = paste0(
                "background-color: #f7f7f7; ",
                "padding: 10px; margin-bottom: 20px;"),
            createNumericInputWithPopover("vaf_het",
                "VAF Heterozygous Threshold:",
                "Define the lower VAF threshold for heterozygous variants.",
                value = 35, min = 0, max = 100, step = 10),
            createNumericInputWithPopover("min_cell_pt",
                "Minimum Cell Percentage:",
                "Enter the minimum percentage of cells with the variant (%).",
                value = 50, min = 0, max = 100),
            createNumericInputWithPopover("min_mut_cell_pt",
                "Minimum Mutated Cell Percentage:",
                "Enter the minimum percentage of mutated cells required (%).",
                value = 1, min = 0, max = 100)),
        actionButton("filter_btn", "Apply Filtering",
            icon = icon("filter"),
            style = paste0(
                "color: #fff; background-color: #337ab7;",
                "border-color: #2e6da4"
            )
        )
    )
}
