"""Small fixtures only: python3 -m unittest discover -s scripts -p 'test_process_species_estimation.py'."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import shlex
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('process', Path(__file__).with_name('process_species_estimation.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ProcessingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.marker = 'S3.1.test'
        self.marker_dir = self.root / self.marker
        self.marker_dir.mkdir()
        self.db = self.root / 'reference database.sdb'
        self.db.mkdir()
        (self.db / 'otus.sqlite3').touch()
        (self.db / 'CONTENTS.json').write_text('{}')
        indices = self.db / 'nucleotide_indices_smafa_naive'
        indices.mkdir()
        (indices / (self.marker + '.smafa_naive_index')).touch()
        self.binary = self.root / 'singlem'
        self.binary.touch()
        self.fasta = self.marker_dir / 'representatives.fasta'
        self.fasta.write_text('>seq1|n_query=1\nAAAA\n'
                              '>otu0|sample2|gene|n_query=0\nCCCC\n'
                              '>seq3|n_query=1\nGGGG\n'
                              '>seq4|n_query=0\nTTTT\n')

    def fake_query(self, cmd, stdout, **kwargs):
        self.assertNotIn('--threads', cmd)
        self.assertIn('--preload-db', cmd)
        self.assertEqual(cmd[cmd.index('--max-divergence') + 1], '2')
        with open(cmd[cmd.index('--query-otu-table') + 1]) as fh:
            queries = list(csv.DictReader(fh, delimiter='\t'))
        stdout.write('query_name\tdivergence\ttaxonomy\n')
        for q in queries:
            if q['sequence'] == 'TTTT':
                continue
            tax = 'd__Bacteria; s__Shared' if q['sequence'] != 'GGGG' else 'd__Bacteria; s__'
            stdout.write(f"{q['sample']}\t2\t{tax}\n")
        return subprocess.CompletedProcess(cmd, 0)

    def run_worker(self):
        m.worker_globdb_marker(str(self.root), self.marker, str(self.db), str(self.binary),
                               str(self.root), chunk_size=1)

    def provenance(self, divergence=2):
        return m.globdb_provenance(str(self.marker_dir), str(self.db), divergence, str(self.binary))

    def test_all_representatives_and_species_collapse_across_batches(self):
        with patch.object(m.subprocess, 'run', side_effect=self.fake_query) as run:
            self.run_worker()
        self.assertEqual(run.call_count, 4)
        summary = json.loads((self.marker_dir / m.GLOBDB_SUMMARY).read_text())
        self.assertEqual(summary['raw_clusters'], 4)
        self.assertEqual(summary['confirmed_clusters'], 3)
        self.assertEqual(summary['confirmed_clusters_with_species'], 2)
        self.assertEqual(summary['distinct_matched_species'], 1)
        self.assertEqual(summary['unmatched_clusters'], 1)
        with open(self.marker_dir / m.GLOBDB_MATCHES) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual(rows[1]['id'], 'otu0|sample2|gene')
        self.assertEqual(rows[2]['globdb_status'], 'confirmed')
        self.assertEqual(rows[2]['species'], '')
        self.assertEqual(rows[3]['globdb_status'], 'no_match')
        self.assertTrue(m.globdb_complete(str(self.marker_dir), self.provenance()))
        self.assertFalse(m.globdb_complete(str(self.marker_dir), self.provenance(1)))
        self.fasta.write_text(self.fasta.read_text() + '>seq5\nACGT\n')
        self.assertFalse(m.globdb_complete(str(self.marker_dir), self.provenance()))

    def test_failed_query_does_not_publish_or_leave_completion(self):
        with patch.object(m.subprocess, 'run', side_effect=self.fake_query):
            self.run_worker()
        original = (self.marker_dir / m.GLOBDB_MATCHES).read_text()
        with patch.object(m.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'singlem')):
            with self.assertRaises(subprocess.CalledProcessError):
                self.run_worker()
        self.assertEqual((self.marker_dir / m.GLOBDB_MATCHES).read_text(), original)
        self.assertFalse(m.globdb_complete(str(self.marker_dir), self.provenance()))

    def test_best_hit_and_first_tie_without_species(self):
        hits = self.root / 'hits.tsv'
        hits.write_text('query_name\tdivergence\ttaxonomy\nq\t2\ts__Far\nq\t1\td__Bacteria\nq\t1\ts__Tie\n')
        self.assertEqual(m._query_best_hits(hits, [('q', 'A')], 2)['q'], (1, m.NA, 'd__Bacteria'))
        with self.assertRaises(ValueError):
            m._query_best_hits(hits, [('other', 'A')], 2)
        hits.write_text('query_name\tdivergence\ttaxonomy\nq\t3\ts__Far\n')
        with self.assertRaises(ValueError):
            m._query_best_hits(hits, [('q', 'A')], 2)

    def test_empty_and_invalid_query_output(self):
        hits = self.root / 'hits.tsv'
        hits.write_text('query_name\tdivergence\ttaxonomy\n')
        self.assertEqual(m._query_best_hits(hits, [('q', 'A')], 2), {})
        hits.write_text('')
        with self.assertRaises(ValueError):
            m._query_best_hits(hits, [('q', 'A')], 2)

    def test_atomic_failure_keeps_existing_output(self):
        path = self.root / 'counts.tsv'
        path.write_text('old')
        with self.assertRaises(RuntimeError):
            with m.atomic_text(path) as fh:
                fh.write('partial')
                raise RuntimeError('failure')
        self.assertEqual(path.read_text(), 'old')

    def test_legacy_histogram_and_correct_new_heading(self):
        path = self.root / 'counts.tsv'
        path.write_text('stratum\tstratum_value\tsample_count\tspecies_count\nyear\t2020\t2\t3\n')
        self.assertEqual(list(m.read_histogram(path)), [('year', '2020', 2, 3)])
        m.write_histogram({('year', '2020', 2): 3}, path)
        self.assertIn('species_count', path.read_text().splitlines()[0])
        self.assertEqual(list(m.read_histogram(path)), [('year', '2020', 2, 3)])

    def test_stratification_deduplicates_same_sample(self):
        (self.root / 'collated_fastas').mkdir()
        (self.root / 'collated_fastas' / (self.marker + '.fasta')).write_text(
            '>otu0|A|gene\nAAAA\n>otu1|A|gene\nCCCC\n>otu0|B|gene\nAAAA\n')
        (self.marker_dir / 'clusters.tsv').write_text('AAAA\tAAAA\nCCCC\tAAAA\n')
        metadata = self.root / 'metadata.tsv'
        metadata.write_text('acc\tyear\torganism\nA\t2020\thuman gut metagenome\nB\t2021\thuman gut metagenome\n')
        m.worker_stratify_marker(str(self.root), self.marker, str(metadata))
        counts = list(m.read_histogram(self.marker_dir / m.STRATIFIED_COUNTS_FILENAME))
        self.assertCountEqual(counts, [
            ('human_all', 'all', 2, 1),
            ('human_per_year', '2020', 1, 1),
            ('human_per_year', '2021', 2, 1),
            ('year', '2020', 1, 1),
            ('year', '2021', 2, 1),
        ])

    def test_year_strata_are_cumulative(self):
        (self.root / 'collated_fastas').mkdir()
        (self.root / 'collated_fastas' / (self.marker + '.fasta')).write_text(
            '>otu0|A|gene\nAAAA\n>otu0|B|gene\nAAAA\n>otu0|C|gene\nCCCC\n>otu1|C|gene\nAAAA\n')
        (self.marker_dir / 'clusters.tsv').write_text('AAAA\tAAAA\nCCCC\tCCCC\n')
        metadata = self.root / 'metadata.tsv'
        metadata.write_text('acc\tyear\torganism\nA\t2018\tsoil\nB\t2020\thuman gut metagenome\n'
                            'C\t2022\tsoil\nD\t2019\tsoil\n')
        m.worker_stratify_marker(str(self.root), self.marker, str(metadata))
        counts = list(m.read_histogram(self.marker_dir / m.STRATIFIED_COUNTS_FILENAME))
        self.assertCountEqual(counts, [
            ('human_all', 'all', 1, 1),
            ('human_per_year', '2020', 1, 1),
            ('human_per_year', '2022', 1, 1),
            ('year', '2018', 1, 1),
            ('year', '2019', 1, 1),
            ('year', '2020', 2, 1),
            ('year', '2022', 3, 1),
            ('year', '2022', 1, 1),
        ])

    def test_submission_requires_ids_and_uses_background(self):
        with patch.object(m.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, stdout='')):
            with self.assertRaises(RuntimeError):
                m._submit_mqsub_batch(['true'], 'test', 1, 1)
        with patch.object(m.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, stdout='qsub stdout: 123.aqua\n')) as run:
            self.assertEqual(m._submit_mqsub_batch(['true'], 'test', 1, 1), ['123.aqua'])
            self.assertIn('--bg', run.call_args.args[0])
            self.assertIn('--no-email', run.call_args.args[0])

    def test_orchestrator_dispatch_summary_and_cached_rerun(self):
        (self.marker_dir / (self.marker + '.tsv')).write_text(
            'id\tn_samples\nseq1\t2\notu0|sample2|gene\t1\nseq3\t1\nseq4\t1\n')
        argv = ['process', str(self.root), '--skip-taxonomy', '--skip-stratify-metadata',
                '--query-db', str(self.db), '--singlem-bin', str(self.binary),
                '--smafa-bin-dir', str(self.root), '--query-chunk-size', '1']
        def submit(cmds, *args, **kwargs):
            for cmd in cmds:
                tokens = shlex.split(cmd)
                with patch.object(sys, 'argv', tokens[1:]):
                    m.main()
            return ['123.aqua'] if cmds else []
        with patch.object(m, '_submit_mqsub_batch', side_effect=submit), \
                patch.object(m, '_wait_for_jobs'), \
                patch.object(m.subprocess, 'run', side_effect=self.fake_query) as query:
            with patch.object(sys, 'argv', argv):
                m.main()
            self.assertEqual(query.call_count, 4)
            self.assertTrue(m.globdb_complete(str(self.marker_dir), self.provenance()),
                            (self.marker_dir / m.GLOBDB_MANIFEST).read_text())
            with patch.object(sys, 'argv', argv):
                m.main()
            self.assertEqual(query.call_count, 4)  # completed queries reused
        with open(self.root / 'sample_count_summary.tsv') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual([(r['sample_count'], r['species_count']) for r in rows], [('1', '3'), ('2', '1')])
        with open(self.root / 'globdb_summary.tsv') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual(rows[0]['distinct_matched_species'], '1')
        self.assertEqual(rows[0]['max_divergence'], '2')

    def test_invalid_batch_size(self):
        with self.assertRaises(ValueError):
            list(m._iter_representative_batches(self.fasta, 0))

    def test_metadata_export_uses_organism_from_sample_attributes(self):
        output = self.root / "metadata.tsv"
        result = subprocess.CompletedProcess([], 0, stdout="", stderr="")
        with patch.object(m.subprocess, "run", return_value=result) as run, \
                patch.object(m.os, "replace"):
            m.export_sample_metadata("database.duckdb", str(output))
        command = run.call_args.args[0]
        query = command[command.index("-c") + 1]
        self.assertIn("m.taxon_name AS organism", query)
        self.assertNotIn("biosample_attributes", query)
        self.assertNotIn("m.organism", query)


if __name__ == '__main__':
    unittest.main()
