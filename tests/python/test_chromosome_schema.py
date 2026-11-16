"""Exercise actual template CSV calls without Nextflow or plotting dependencies.

Run: python -m unittest discover -s tests/python -v
Requires the module environment's polars (2.0.0) and pandas.
"""

import ast
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

import pandas as pd
import polars as pl


ROOT = Path(__file__).resolve().parents[2]
BED = ["chr", "start", "end", "name", "score", "strand"]
META_BED = BED + ["sample", "tool"]
ANNOTATION = BED + [
    "thickStart", "thickEnd", "itemRgb", "exonCount", "exonSizes",
    "exonOffsets", "readNumber", "circType", "geneName", "isoformName",
    "index", "flankIntron",
]

# Template, call index, file columns, chromosome/start/end column names.
CASES = [
    ("combinebeds/filter/templates/filter.py", 0, META_BED, "chr", "start", "end"),
    ("combinebeds/counts/templates/counts.py", 0, BED, "chr", "start", "end"),
    ("combinebeds/counts/templates/counts.py", 1, META_BED, "chr", "start", "end"),
    ("combinebeds/reads/templates/reads.py", 0, BED + ["reads"], "chr", "start", "end"),
    ("combinebeds/shifts/templates/shifts.py", 0, META_BED, "chr", "start", "end"),
    ("jccirc/prep/templates/prep.py", 0, ANNOTATION, "chr", "start", "end"),
    ("jccirc/prep/templates/prep.py", 1, BED + ["reads"], "chr", "start", "end"),
    ("ciri/build_list/templates/build_list.py", 0, [f"column_{i}" for i in range(12)], "column_0", "column_1", "column_2"),
    ("ciri/build_list/templates/build_list.py", 1, ["BSJ", "Chr", "Start", "End"], "Chr", "Start", "End"),
    ("circtools/unify/templates/unify.py", 0, ["Chr", "Start", "End", "paired.junctions"], "Chr", "Start", "End"),
    ("circtools/unify/templates/unify.py", 1, ["Chr", "Start", "End", "Strand"], "Chr", "Start", "End"),
    ("ciri/cirivis_unify/templates/cirivis_unify.py", 0, ["Chr", "start", "end"], "Chr", "start", "end"),
]


def reader_call(relative, index):
    source = (ROOT / "modules/local" / relative).read_text()
    # Nextflow interpolates placeholders and unescapes doubled backslashes.
    source = re.sub(r"\$\{[^}]*\}", "fixture", source).replace("\\\\", "\\")
    tree = ast.parse(source)
    calls = sorted(
        (node for node in ast.walk(tree) if isinstance(node, ast.Call)
         and isinstance(node.func, ast.Attribute)
         and node.func.attr in {"read_csv", "scan_csv"}),
        key=lambda node: node.lineno,
    )
    return calls[index]


class ChromosomeSchemaTests(unittest.TestCase):
    def test_coordinate_readers(self):
        for relative, index, columns, chrom, start, end in CASES:
            # Numeric-only files must also retain chromosome text for joins.
            for chromosomes in (["1"] * 101, ["1"] * 101 + ["X"]):
                with self.subTest(template=relative, call=index, last=chromosomes[-1]):
                    call = reader_call(relative, index)
                    context = {"pl": pl, "pd": pd, "columns": BED}
                    if "combinebeds/reads/" in relative:
                        context["columns"] = BED + ["reads"]
                    keywords = {
                        kw.arg: eval(compile(ast.Expression(kw.value), relative, "eval"), context)
                        for kw in call.keywords
                    }
                    has_header = keywords.get("has_header", True) if call.func.value.id == "pl" else True
                    with tempfile.TemporaryDirectory() as directory:
                        path = Path(directory) / "input.bed"
                        rows = ["\t".join(columns)] if has_header else []
                        for chromosome in chromosomes:
                            values = {chrom: chromosome, start: "100", end: "200"}
                            rows.append("\t".join(values.get(col, ".") for col in columns))
                        path.write_text("\n".join(rows) + "\n")
                        reader = getattr(context[call.func.value.id], call.func.attr)
                        frame = reader(str(path), **keywords)
                        if isinstance(frame, pl.LazyFrame):
                            frame = frame.collect()
                        if isinstance(frame, pl.DataFrame):
                            self.assertEqual(frame.schema[chrom], pl.String)
                            self.assertTrue(frame.schema[start].is_integer())
                            self.assertTrue(frame.schema[end].is_integer())
                            self.assertEqual(frame[chrom].to_list(), chromosomes)
                            self.assertEqual(frame[start].to_list(), [100] * len(chromosomes))
                            self.assertEqual(frame[end].to_list(), [200] * len(chromosomes))
                        else:
                            self.assertTrue(pd.api.types.is_string_dtype(frame[chrom]))
                            self.assertTrue(pd.api.types.is_integer_dtype(frame[start]))
                            self.assertTrue(pd.api.types.is_integer_dtype(frame[end]))
                            self.assertEqual(frame[chrom].tolist(), chromosomes)
                            self.assertEqual(frame[start].tolist(), [100] * len(chromosomes))
                            self.assertEqual(frame[end].tolist(), [200] * len(chromosomes))

    def test_ciri_consensus_annotation_join(self):
        # Exercise the entire template, including BED12 projection, coordinate
        # conversion and annotation matching, beyond the isolated reader seam.
        template = ROOT / "modules/local/ciri/build_list/templates/build_list.py"
        chromosomes = [str(i) for i in range(1, 102)] + ["X"]
        for max_shift in (0, 1):
            with self.subTest(max_shift=max_shift), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                bed = directory / "consensus.bed"
                anno = directory / "annotation.tsv"
                bed.write_text("".join(
                    f"{chrom}\t100\t200\tcirc{chrom}\t5\t+\t100\t200\t0\t1\t100\t0\n"
                    for chrom in chromosomes
                ))
                anno.write_text("BSJ\tChr\tStart\tEnd\n" + "".join(
                    f"circ{chrom}\t{chrom}\t{101 + max_shift}\t200\n" for chrom in chromosomes
                ))
                values = {
                    "consensus": str(bed), "ciri_annotation": str(anno),
                    "max_shift": str(max_shift), "prefix": "result", "task.process": "BUILD_LIST",
                }
                source = re.sub(r"\$\{([^}]*)\}", lambda match: values[match[1]], template.read_text())
                source = source.replace("\\\\", "\\")
                script = directory / "build_list.py"
                script.write_text(source)
                result = subprocess.run([sys.executable, str(script)], cwd=directory, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(
                    sorted((directory / "result.txt").read_text().splitlines()),
                    sorted(f"circ{chrom}" for chrom in chromosomes),
                )


if __name__ == "__main__":
    unittest.main()
