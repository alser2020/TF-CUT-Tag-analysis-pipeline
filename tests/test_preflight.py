from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import yaml

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "preflight.py"
SPEC = importlib.util.spec_from_file_location("cuttag_preflight", MODULE_PATH)
preflight = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(preflight)


class PreflightTests(unittest.TestCase):
    def make_project(self, root: Path, *, bad_gtf=False, duplicate_sample=False):
        (root / "data").mkdir()
        for sample in ("TF_R1", "TF_R2", "IgG_R1", "IgG_R2"):
            (root / "data" / f"{sample}_1.fq.gz").write_bytes(b"")
            (root / "data" / f"{sample}_2.fq.gz").write_bytes(b"")
        (root / "genome.fa").write_text(">1\nACGTACGT\n", encoding="utf-8")
        gtf_contig = "2" if bad_gtf else "1"
        (root / "genes.gtf").write_text(
            f'{gtf_contig}\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "g1";\n',
            encoding="utf-8",
        )
        rows = [
            ("TF_R1", "TF", 1),
            ("TF_R2", "TF", 2),
            ("IgG_R1", "IgG", 1),
            ("IgG_R2", "IgG", 2),
        ]
        if duplicate_sample:
            rows[-1] = ("TF_R1", "IgG", 2)
        with (root / "samples.tsv").open("w", encoding="utf-8") as handle:
            handle.write("sample_id\tgroup_id\treplicate\tlabel\tr1\tr2\n")
            for sample, group, replicate in rows:
                handle.write(
                    f"{sample}\t{group}\t{replicate}\t{sample}\t"
                    f"data/{sample}_1.fq.gz\tdata/{sample}_2.fq.gz\n"
                )
        (root / "comparisons.tsv").write_text(
            "comparison_id\ttreatment_group\tcontrol_group\nTF_vs_IgG\tTF\tIgG\n",
            encoding="utf-8",
        )
        config = {
            "project": {
                "name": "test",
                "output_dir": "results",
                "samples": "samples.tsv",
                "comparisons": "comparisons.tsv",
            },
            "reference": {
                "assembly": "test_assembly",
                "fasta": "genome.fa",
                "gtf": "genes.gtf",
                "bowtie2_index": None,
                "effective_genome_size": 8,
                "chromosomes": ["1"],
            },
            "analysis": {
                "fastp": {
                    "qualified_quality_phred": 15,
                    "length_required": 36,
                    "cut_front": False,
                    "cut_tail": True,
                    "cut_window_size": 4,
                    "cut_mean_quality": 20,
                    "average_qual": 0,
                },
                "alignment": {
                    "mode": "end-to-end",
                    "preset": "very-sensitive",
                    "no_mixed": True,
                    "no_discordant": True,
                    "insert_min": 10,
                    "insert_max": 700,
                    "mapq_min": 30,
                },
                "coverage": {
                    "bin_size": 25,
                    "smooth_length": 50,
                    "extend_reads": 200,
                    "normalize_using": "RPGC",
                    "ignore_duplicates": False,
                },
                "qc": {
                    "summary_bin_size": 25,
                    "correlation_method": "spearman",
                    "fingerprint_bin_size": 25,
                    "fingerprint_samples": 1000,
                    "tss_before": 3000,
                    "tss_after": 3000,
                    "tss_bin_size": 25,
                },
                "macs3": {
                    "format": "BAMPE",
                    "qvalue": 0.05,
                    "broad": False,
                    "call_summits": True,
                    "trackline": False,
                    "keep_dup": "auto",
                },
                "idr": {"mode": "auto", "threshold": 0.05, "rank": "p.value"},
                "annotation": {
                    "enabled": False,
                    "method": "homer",
                    "promoter_max_abs_tss_distance": 5000,
                },
                "motif": {
                    "enabled": True,
                    "methods": ["homer", "memechip"],
                    "input_scope": "all",
                    "homer": {"size": 200, "lengths": [6, 8, 10, 12]},
                    "memechip": {"flank": 100, "center_cut": 0, "min_width": 6, "max_width": 12, "meme_p": 16, "meme_nmotifs": 3, "streme_nmotifs": 0, "known_motif_db": None},
                },
            },
            "outputs": {"keep_trimmed_reads": True},
        }
        config_path = root / "config.yaml"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        return config_path

    def test_release_bin_and_center_cut_defaults(self):
        release = Path(__file__).resolve().parents[1]
        config = yaml.safe_load((release / "config/config.example.yaml").read_text())
        analysis = config["analysis"]
        self.assertEqual(analysis["coverage"]["bin_size"], 25)
        for key in ("summary_bin_size", "fingerprint_bin_size", "tss_bin_size"):
            self.assertEqual(analysis["qc"][key], 25)
        self.assertEqual(analysis["motif"]["memechip"]["center_cut"], 0)
        self.assertEqual(analysis["motif"]["memechip"]["streme_nmotifs"], 0)

    def test_valid_project(self):
        with tempfile.TemporaryDirectory() as tmp:
            config_path = self.make_project(Path(tmp))
            report, normalized = preflight.validate(config_path)
            self.assertEqual(report.errors, [])
            self.assertEqual(normalized["sample_count"], 4)
            self.assertTrue(any("将运行 IDR" in message for message in report.info))

    def test_reference_mismatch_is_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            config_path = self.make_project(Path(tmp), bad_gtf=True)
            report, _ = preflight.validate(config_path)
            self.assertTrue(any("GTF" in message and "不在 FASTA" in message for message in report.errors))

    def test_duplicate_sample_is_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            config_path = self.make_project(Path(tmp), duplicate_sample=True)
            report, _ = preflight.validate(config_path)
            self.assertTrue(any("sample_id 重复" in message for message in report.errors))

    def test_three_replicates_warn_and_skip_idr(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = self.make_project(root)
            (root / "data" / "TF_R3_1.fq.gz").write_bytes(b"")
            (root / "data" / "TF_R3_2.fq.gz").write_bytes(b"")
            with (root / "samples.tsv").open("a", encoding="utf-8") as handle:
                handle.write("TF_R3\tTF\t3\tTF_R3\tdata/TF_R3_1.fq.gz\tdata/TF_R3_2.fq.gz\n")
            report, _ = preflight.validate(config_path)
            self.assertEqual(report.errors, [])
            self.assertTrue(any("有 3 个重复" in message and "跳过 IDR" in message for message in report.warnings))

    def test_promoter_motif_requires_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = self.make_project(root)
            config = yaml.safe_load(config_path.read_text())
            config["analysis"]["motif"]["input_scope"] = "promoter"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False))
            report, _ = preflight.validate(config_path)
            self.assertTrue(any("必须启用 annotation" in message for message in report.errors))

    def test_invalid_qvalue_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = self.make_project(root)
            config = yaml.safe_load(config_path.read_text())
            config["analysis"]["macs3"]["qvalue"] = 2
            config_path.write_text(yaml.safe_dump(config, sort_keys=False))
            report, _ = preflight.validate(config_path)
            self.assertTrue(any("qvalue" in message for message in report.errors))


if __name__ == "__main__":
    unittest.main()
