from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path


class FilterHomerAnnotationTests(unittest.TestCase):
    def test_distance_filter(self):
        root = Path(__file__).resolve().parents[1]
        script = root / "scripts" / "filter_homer_annotation.py"
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            annotation = tmp / "annotation.tsv"
            summits = tmp / "summits.bed"
            bed = tmp / "promoter.bed"
            selected = tmp / "promoter_summits.bed"
            annotation.write_text(
                "PeakID\tChr\tStart\tEnd\tStrand\tPeak Score\tFocus Ratio\tAnnotation\tDistance to TSS\n"
                "peak1\t1\t10\t11\t+\t10\t0\tIntergenic\t4000\n"
                "peak2\t1\t20\t21\t+\t10\t0\tIntergenic\t9000\n",
                encoding="utf-8",
            )
            summits.write_text(
                "1\t10\t11\tpeak1\t10\n1\t20\t21\tpeak2\t10\n",
                encoding="utf-8",
            )
            subprocess.run(
                [
                    "python3", str(script),
                    "--annotation", str(annotation),
                    "--summits", str(summits),
                    "--max-distance", "5000",
                    "--bed", str(bed),
                    "--summit-bed", str(selected),
                ],
                check=True,
            )
            self.assertIn("peak1", bed.read_text())
            self.assertNotIn("peak2", bed.read_text())
            self.assertIn("peak1", selected.read_text())
            self.assertNotIn("peak2", selected.read_text())


if __name__ == "__main__":
    unittest.main()
