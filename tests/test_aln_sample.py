import pytest
from sat.scripts.aln_sample import aln_sample_main


class Args:
    pass


class TestAlnSampleMain:

    def test_random_sample_with_header(self, tmp_path):
        """Test random sampling using the real test data file."""
        args = Args()
        args.alignment_file = (
            "tests/test_data/foldseek_related/top_query_per_cluster_tax.m8"
        )
        args.output_file = f"{tmp_path}/sampled.m8"
        args.alignment_fields = ""
        args.n_alignments = 2
        args.query_column = "query"
        args.sort_column = ""
        args.random_seed = 42

        aln_sample_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # First line should be header
        assert lines[0].startswith("query\t")

        # Count alignments per query
        header = lines[0].strip().split("\t")
        query_idx = header.index("query")
        query_counts = {}
        for line in lines[1:]:
            parts = line.strip().split("\t")
            q = parts[query_idx]
            query_counts[q] = query_counts.get(q, 0) + 1

        # No query should have more than 2 alignments
        for q, count in query_counts.items():
            assert count <= 2, f"Query {q} has {count} alignments, expected <= 2"

    def test_random_sample_reproducibility(self, tmp_path):
        """Same seed produces identical output."""
        args1 = Args()
        args1.alignment_file = (
            "tests/test_data/foldseek_related/top_query_per_cluster_tax.m8"
        )
        args1.output_file = f"{tmp_path}/sampled1.m8"
        args1.alignment_fields = ""
        args1.n_alignments = 2
        args1.query_column = "query"
        args1.sort_column = ""
        args1.random_seed = 42

        args2 = Args()
        args2.alignment_file = args1.alignment_file
        args2.output_file = f"{tmp_path}/sampled2.m8"
        args2.alignment_fields = ""
        args2.n_alignments = 2
        args2.query_column = "query"
        args2.sort_column = ""
        args2.random_seed = 42

        aln_sample_main(args1)
        aln_sample_main(args2)

        with open(args1.output_file) as f:
            content1 = f.read()
        with open(args2.output_file) as f:
            content2 = f.read()

        assert content1 == content2

    def test_top_n_by_sort_column(self, tmp_path):
        """Top-N selection picks the highest values."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tscore\n"
            "q1\tt1\t0.9\n"
            "q1\tt2\t0.5\n"
            "q1\tt3\t0.7\n"
            "q2\tt4\t0.3\n"
            "q2\tt5\t0.8\n"
        )

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.n_alignments = 1
        args.query_column = "query"
        args.sort_column = "score"
        args.random_seed = 42

        aln_sample_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Header + 2 rows (one per query)
        assert len(lines) == 3
        assert lines[0].strip() == "query\ttarget\tscore"

        data = [line.strip().split("\t") for line in lines[1:]]
        data_dict = {row[0]: row for row in data}

        # q1 should keep t1 (score 0.9), q2 should keep t5 (score 0.8)
        assert data_dict["q1"][1] == "t1"
        assert data_dict["q2"][1] == "t5"

    def test_fewer_than_n_keeps_all(self, tmp_path):
        """When N exceeds group size, all rows are kept."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tscore\n"
            "q1\tt1\t0.9\n"
            "q1\tt2\t0.5\n"
            "q2\tt3\t0.3\n"
        )

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.n_alignments = 10
        args.query_column = "query"
        args.sort_column = ""
        args.random_seed = 42

        aln_sample_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Header + all 3 data rows
        assert len(lines) == 4

    def test_explicit_alignment_fields(self, tmp_path):
        """Headerless input with provided column names works."""
        input_file = tmp_path / "input_no_header.m8"
        input_file.write_text(
            "q1\tt1\t0.9\n"
            "q1\tt2\t0.5\n"
            "q2\tt3\t0.3\n"
        )

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = "query,target,score"
        args.n_alignments = 1
        args.query_column = "query"
        args.sort_column = "score"
        args.random_seed = 42

        aln_sample_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Header + 2 rows (one per query)
        assert len(lines) == 3
        assert lines[0].strip() == "query\ttarget\tscore"

    def test_custom_query_column(self, tmp_path):
        """Non-default query column name works."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "protein_id\ttarget\tscore\n"
            "p1\tt1\t0.9\n"
            "p1\tt2\t0.5\n"
            "p1\tt3\t0.7\n"
            "p2\tt4\t0.3\n"
        )

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.n_alignments = 2
        args.query_column = "protein_id"
        args.sort_column = ""
        args.random_seed = 42

        aln_sample_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # p1 has 3 rows, should be sampled to 2. p2 has 1, kept as-is.
        assert len(lines) == 4  # header + 3 data rows

    def test_missing_query_column_raises(self, tmp_path):
        """ValueError when query_column doesn't exist in the data."""
        input_file = tmp_path / "input.m8"
        input_file.write_text("query\ttarget\tscore\nq1\tt1\t0.5\n")

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.n_alignments = 1
        args.query_column = "nonexistent"
        args.sort_column = ""
        args.random_seed = 42

        with pytest.raises(ValueError, match="Query column 'nonexistent' not found"):
            aln_sample_main(args)

    def test_missing_sort_column_raises(self, tmp_path):
        """ValueError when sort_column doesn't exist in the data."""
        input_file = tmp_path / "input.m8"
        input_file.write_text("query\ttarget\tscore\nq1\tt1\t0.5\n")

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.n_alignments = 1
        args.query_column = "query"
        args.sort_column = "nonexistent"
        args.random_seed = 42

        with pytest.raises(ValueError, match="Sort column 'nonexistent' not found"):
            aln_sample_main(args)

    def test_scientific_notation_sort(self, tmp_path):
        """Scientific notation values are handled correctly in sort."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tscore\n"
            "q1\tt1\t4.591E-01\n"
            "q1\tt2\t9.510E-04\n"
            "q1\tt3\t6.000E-01\n"
        )

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.n_alignments = 2
        args.query_column = "query"
        args.sort_column = "score"
        args.random_seed = 42

        aln_sample_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Should keep 2 highest: 6.000E-01 and 4.591E-01
        assert len(lines) == 3
        targets = [line.strip().split("\t")[1] for line in lines[1:]]
        assert "t3" in targets  # 0.6
        assert "t1" in targets  # ~0.459
        assert "t2" not in targets  # ~0.00095 dropped
