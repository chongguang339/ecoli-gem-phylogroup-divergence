from pathlib import Path

from docx import Document


DOCX = Path(r"D:\pipeline\analysis_v5\02.论文\Evolutionary Divergence of Genomic Defense.docx")


def set_text(paragraph, text):
    paragraph.clear()
    paragraph.add_run(text)


def main():
    doc = Document(DOCX)
    for para in doc.paragraphs:
        text = para.text.strip()
        if text.startswith("The four available curated-reference genomes"):
            set_text(
                para,
                "The reference-model structural benchmark comprised five model systems (Figure 1D). Four systems (iB21_1397, iEC1364_W, iECDH10B_1368, and iJN1463) had independently logged paired default/optimized runtime records and showed substantial acceleration, with a median speedup of 12.76-fold and a range of 11.30-15.52-fold. MG1655/iML1515 was retained as the fifth structural reference system and used the benchmark-average runtime basis for five-model summaries. Across the four independently timed reference systems, optimized reaction retention ranged from 77.0% to 87.5%. When reactions present in the default reconstruction but absent after optimization were compared against the matched BiGG reference models, only 10.55-17.53% of these removed reactions were present in the curated references. Thus, most removed reactions were not supported by the corresponding reference reconstructions, although the reference-supported subset of removed reactions should be audited before making strain-specific functional claims.",
            )
        elif text.startswith("Biomass flux values among viable models"):
            set_text(
                para,
                "Biomass flux values among viable models range from 0.526 to 0.584 mmol gDW-1 h-1 under standardized rich medium conditions, with median flux = 0.5826 mmol gDW-1 h-1, closely matching the MG1655 reference value (0.5826 mmol gDW-1 h-1). This consistency indicates that the optimized workflow produced solvable biomass equations under the standardized rich-medium test condition, but it should not be treated as proof of broad substrate-level predictive performance; substrate-specific behavior is better assessed by the carbon-source benchmarks described above.",
            )
        elif text.startswith("The computational benchmark uses a distribution-based assessment"):
            set_text(
                para,
                "The computational benchmark uses a distribution-based assessment across 20 N50-stratified low-contiguity assemblies annotated with phylogroup and CDS metadata (A (n = 4), B1 (n = 7), B2 (n = 4), C (n = 1), E (n = 2), cryptic/U (n = 2); CDS_Count range 4,252-5,584). Optimized gapseq reduced median runtime from 19,821.5 s to 1,260.5 s, corresponding to a median 15.81-fold speedup and a narrow IQR of 15.40-16.35. Five reference model systems were organized for structure and model-consistency benchmarking; MG1655/iML1515 used the same benchmark-average runtime basis, giving a five-model reference interpretation centered on the 12.76-fold speedup estimate while keeping measured 20-genome runtime distributions separate. This distributional benchmark is useful for cohort-scale planning but should not be treated as direct evidence of universal substrate-level predictive performance.",
            )
        elif text.startswith("In conclusion, this study presents"):
            set_text(
                para,
                "In conclusion, this study presents an optimized and systematically benchmarked E. coli GEM reconstruction framework applied to a 1,000-genome phylogroup-matched cohort. Stratified benchmarking across 20 assemblies showed a median 15.81-fold runtime acceleration, paired XML reaction matrices identified 111.0 optimized-only reactions per genome at the median, and five reference model systems supported reaction-set and model-consistency benchmarking, with MG1655/iML1515 using the same benchmark-average runtime basis. Automated SD01 carbon-source mapping translated 738/770 aerobic/anaerobic environments into model-readable exchanges and produced reproducible FBA prediction matrices for curated BiGG and gapseq models, providing a transparent in silico consistency benchmark rather than independent experimental phenotype validation. FastANI residual analysis showed that KEGG pathway composition is strongly tied to continuous core-genome distance, whereas discrete phylogroup labels retain smaller residual effects for KEGG, COG, GEM topology, CRISPR spacer composition, CRISPR target burden, and SignalP/TMHMM surface architecture. The main validation priority is true multi-strain experimental Biolog or carbon-source data for substrate-level predictive benchmarking.",
            )
    doc.save(DOCX)
    print(f"saved {DOCX}")


if __name__ == "__main__":
    main()
