import xml.etree.ElementTree as ET
from pathlib import Path


FILES = [
    Path(r"D:\pipeline\analysis_v4\data\model\gapseq默认参数模型\GCA_000005845.2_ASM584v2_genomic.xml"),
    Path(r"D:\pipeline\analysis_v4\data\model\diamond_enterobacterales_fastdb\GCF_000005845.2_ASM584v2_genomic.xml"),
    Path(r"D:\pipeline\analysis_v4\data\model\iML1515.xml"),
]


def reaction_ids(path):
    root = ET.parse(path).getroot()
    ids = set()
    for elem in root.iter():
        if elem.tag.endswith("reaction"):
            rid = elem.attrib.get("id")
            if rid:
                ids.add(rid)
    return ids


def main():
    sets = {}
    for path in FILES:
        ids = reaction_ids(path)
        sets[path.name] = ids
        print(path.name, len(ids))
    default = sets[FILES[0].name]
    optimized = sets[FILES[1].name]
    curated = sets[FILES[2].name]
    print("optimized_only", len(optimized - default))
    print("default_only", len(default - optimized))
    print("shared", len(default & optimized))
    print("jaccard", len(default & optimized) / len(default | optimized))
    print("optimized_retention", len(optimized) / len(default))
    print("default_removed_supported_by_iML1515", len((default - optimized) & curated))
    print("optimized_only_supported_by_iML1515", len((optimized - default) & curated))


if __name__ == "__main__":
    main()
