from __future__ import annotations

from pathlib import Path
import re

import pandas as pd
from docx import Document


ROOT = Path(r"D:\pipeline\analysis_v5")
DOCX = ROOT / "02.论文" / "Evolutionary Divergence of Genomic Defense.docx"
OUT = ROOT / "03.analysis" / "manuscript_consistency_audit.md"
CSV = ROOT / "03.analysis" / "manuscript_consistency_audit_findings.csv"


def ptxt(doc: Document, needle: str) -> list[tuple[int, str]]:
    out = []
    for i, p in enumerate(doc.paragraphs):
        if needle.lower() in p.text.lower():
            out.append((i, p.text))
    return out


def add_finding(rows, severity, topic, finding, evidence, recommendation):
    rows.append(
        {
            "Severity": severity,
            "Topic": topic,
            "Finding": finding,
            "Evidence": evidence,
            "Recommendation": recommendation,
        }
    )


def main() -> None:
    doc = Document(DOCX)
    rows = []
    lines = ["# Manuscript Consistency Audit", ""]

    # Cohort metrics
    ann = pd.read_csv(ROOT / "00.data" / "recalculated_annotations.csv")
    cog = pd.read_csv(ROOT / "00.data" / "genome_assembly_and_cog_metrics.csv")
    cohort_counts = ann["Phylogroup"].value_counts().to_dict()
    expected = {"A": 298, "B2": 244, "B1": 200, "D": 103, "E": 73, "C": 36, "F": 31, "G": 11, "U": 4}
    if cohort_counts != expected:
        add_finding(
            rows,
            "Major",
            "Cohort counts",
            "Phylogroup counts in recalculated_annotations.csv differ from the stated revised cohort.",
            f"Observed {cohort_counts}; expected {expected}.",
            "Reconcile `recalculated_annotations.csv` and manuscript counts before submission.",
        )
    else:
        lines.append(f"- Cohort counts match revised design: {expected}.")

    if set(ann["Genome_ID"]) != set(cog["Genome_ID"]):
        add_finding(
            rows,
            "Major",
            "Genome linkage",
            "Genome_ID sets differ between recalculated_annotations.csv and genome_assembly_and_cog_metrics.csv.",
            f"annotations={ann['Genome_ID'].nunique()}, COG={cog['Genome_ID'].nunique()}, overlap={len(set(ann['Genome_ID']) & set(cog['Genome_ID']))}",
            "Use a single joined cohort table for all downstream analyses or document exclusions.",
        )
    else:
        lines.append("- Genome_ID sets match between annotation and assembly/COG matrices.")

    # Benchmark metrics
    bm = pd.read_csv(ROOT / "03.analysis" / "benchmark_20_with_phylogroup_cds_reaction_delta.csv")
    bm_counts = bm["Phylogroup"].value_counts().to_dict()
    lines.append(f"- Benchmark metadata table has {len(bm)} genomes; phylogroups: {bm_counts}; CDS range {bm['CDS_Count'].min()}-{bm['CDS_Count'].max()}.")
    if {"D", "F", "G"} & set(bm_counts):
        add_finding(
            rows,
            "Minor",
            "Benchmark scope wording",
            "Benchmark table includes D/F/G despite manuscript scope language saying they are absent.",
            str(bm_counts),
            "Update benchmark wording.",
        )
    if not {"D", "F", "G"} & set(bm_counts):
        lines.append("- Benchmark scope language is appropriate: D/F/G are absent from the 20-genome stress-test subset.")

    # CRISPR target metrics
    target_summary = pd.read_csv(ROOT / "03.analysis" / "crispr_target_matching_analysis" / "crispr_target_match_summary.csv")
    all_eff = pd.read_csv(ROOT / "03.analysis" / "crispr_target_matching_analysis" / "crispr_target_matching_phylogroup_effect_sizes_all_genomes.csv")
    primary = target_summary[target_summary["Filter"] == "mismatch_le2"].set_index("Target_Type")
    lines.append(
        "- CRISPR target matching primary filter yields "
        f"{int(primary.loc['phage','Unique_targeted_spacers'])} phage-targeting and "
        f"{int(primary.loc['plasmid','Unique_targeted_spacers'])} plasmid-targeting unique spacers."
    )
    any_eff = all_eff.set_index("Metric").loc["Any_targeted_occurrences"]
    plasmid_eff = all_eff.set_index("Metric").loc["Plasmid_targeted_occurrences"]
    lines.append(f"- CRISPR all-genome target burden effect sizes: plasmid eta2={plasmid_eff.eta2:.3f}, any-target eta2={any_eff.eta2:.3f}.")

    # KEGG/PERMANOVA metrics
    kegg = pd.read_csv(ROOT / "03.analysis" / "pathway_hellinger_permanova_full.csv")
    permdisp = pd.read_csv(ROOT / "03.analysis" / "pathway_hellinger_permdisp_full.csv")
    kegg_r2 = float(kegg["R2"].iloc[0]) if "R2" in kegg.columns else None
    kegg_p_col = "p_value" if "p_value" in kegg.columns else "p"
    disp_p_col = "p_value" if "p_value" in permdisp.columns else "p"
    kegg_p = float(kegg[kegg_p_col].iloc[0]) if kegg_p_col in kegg.columns else None
    disp_p = float(permdisp[disp_p_col].iloc[0]) if disp_p_col in permdisp.columns else None
    lines.append(f"- KEGG Hellinger PERMANOVA file reports R2={kegg_r2}, p={kegg_p}; PERMDISP p={disp_p}.")

    # ANI residuals
    ani = pd.read_csv(ROOT / "03.analysis" / "ani_residual_phylogroup_effects.csv")
    if "KEGG" in " ".join(ani.astype(str).agg(" ".join, axis=1)):
        lines.append("- ANI residual phylogroup effects file is present and includes functional matrices.")
    else:
        add_finding(
            rows,
            "Major",
            "ANI residual analysis",
            "ANI residual output could not be clearly matched to named functional matrices.",
            f"Columns: {list(ani.columns)}",
            "Check `ani_residual_phylogroup_effects.csv` labels and manuscript wording.",
        )

    # Text-level risks
    text = "\n".join(p.text for p in doc.paragraphs)
    risky = ["explains 56.1", "R² = 0.561", "R2 = 0.561", "B2 (n = 1)", "A (n = 405)", "no loss of final reaction content"]
    found_risky = [r for r in risky if r in text]
    if found_risky:
        add_finding(
            rows,
            "Critical",
            "Obsolete claims",
            "Old high-risk claims remain in the manuscript.",
            "; ".join(found_risky),
            "Remove or replace obsolete statements.",
        )
    else:
        lines.append("- Old high-risk claims were not found in the current DOCX text.")

    # Known missing publication metadata
    if "[repository URL]" in text or "[Funding details]" in text:
        add_finding(
            rows,
            "Major",
            "Submission metadata",
            "Placeholder repository/funding text remains in declarations.",
            "Found `[repository URL]` and/or `[Funding details]` in DOCX.",
            "Replace with final repository URL/DOI and funding statement before submission.",
        )

    # Supplementary figure citation consistency
    if "Supplementary Figure S1" in text and not (ROOT / "01.figures" / "FigureS1_crispr_target_matching.png").exists():
        add_finding(
            rows,
            "Major",
            "Supplementary figure",
            "Manuscript cites Supplementary Figure S1, but the figure file is missing.",
            "Expected 01.figures/FigureS1_crispr_target_matching.png.",
            "Generate or remove the citation.",
        )
    elif "Supplementary Figure S1" in text:
        lines.append("- Supplementary Figure S1 is cited and the PNG exists.")

    # Subtle logic risks
    best = pd.read_csv(ROOT / "03.analysis" / "crispr_target_matching_analysis" / "crispr_target_hits_best_mismatch_le2.csv")
    phage_ids = set(best.loc[best["Target_Type"] == "phage", "Unique_Spacer_ID"])
    plasmid_ids = set(best.loc[best["Target_Type"] == "plasmid", "Unique_Spacer_ID"])
    overlap_n = len(phage_ids & plasmid_ids)
    if "overlapped for 39 unique spacers" in text and "database-match categories" in text:
        lines.append(f"- CRISPR phage/plasmid target overlap is acknowledged in text: {overlap_n} shared unique spacers.")
    else:
        add_finding(
            rows,
            "Minor",
            "CRISPR target interpretation",
            "Phage-targeted and plasmid-targeted spacers are not necessarily mutually exclusive because mobile-element databases can share homologous sequence.",
            f"Observed overlap: {overlap_n} Unique_Spacer_IDs.",
            "Add a sentence in Methods/Results or Supplementary notes that target categories are database-match categories and may overlap.",
        )
    if "In this dataset, phylogroup C" in text and "n = 36" in text:
        lines.append("- Small-group CRISPR target interpretation is worded cautiously for phylogroup C.")
    else:
        add_finding(
            rows,
            "Minor",
            "Small phylogroups",
            "C, G, and U subgroup estimates are based on small n, and C shows 100% plasmid-targeted genomes.",
            "C n=36, G n=11, U n=4; extreme fractions are plausible but sensitive to sample size.",
            "Use cautious wording such as 'in this dataset' and avoid broad biological generalization for small groups.",
        )
    add_finding(
        rows,
        "Major",
        "Phenotype validation",
        "The main remaining validation gap is still true multi-strain experimental Biolog/carbon-source benchmarking.",
        "Manuscript acknowledges MG1655 proof-of-concept and lack of independent multi-strain truth data.",
        "Do not claim publication-ready functional accuracy beyond single-strain proof-of-concept unless experimental truth data are added.",
    )

    findings = pd.DataFrame(rows)
    findings.to_csv(CSV, index=False, encoding="utf-8-sig")

    lines.extend(["", "## Findings", ""])
    if findings.empty:
        lines.append("No inconsistencies found.")
    else:
        for row in findings.itertuples(index=False):
            lines.append(f"### {row.Severity}: {row.Topic}")
            lines.append(f"- Finding: {row.Finding}")
            lines.append(f"- Evidence: {row.Evidence}")
            lines.append(f"- Recommendation: {row.Recommendation}")
            lines.append("")

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT}")
    print(f"Wrote {CSV}")
    print(findings[["Severity", "Topic", "Finding"]].to_string(index=False))


if __name__ == "__main__":
    main()
