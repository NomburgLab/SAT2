import pytest
from sat.scripts.aln_filter import aln_filter_main, parse_filter_value


class TestParseFilterValue:
    """Tests for the parse_filter_value helper function."""

    def test_parse_regular_float(self):
        assert parse_filter_value("0.5") == 0.5
        assert parse_filter_value("1.0") == 1.0

    def test_parse_scientific_notation_lowercase(self):
        assert parse_filter_value("4.591e-01") == pytest.approx(0.4591)
        assert parse_filter_value("1.2e-12") == pytest.approx(1.2e-12)

    def test_parse_scientific_notation_uppercase(self):
        assert parse_filter_value("4.591E-01") == pytest.approx(0.4591)
        assert parse_filter_value("9.510E-04") == pytest.approx(0.000951)

    def test_parse_invalid_value(self):
        with pytest.raises(ValueError, match="Cannot convert"):
            parse_filter_value("invalid")


class TestAlnFilterMain:
    """Tests for the main aln_filter_main function."""

    def test_basic_filtering_with_header(self, tmp_path):
        """Test basic filtering when input file has a header."""

        class Args:
            pass

        args = Args()
        args.alignment_file = (
            "tests/test_data/foldseek_related/top_query_per_cluster_tax.m8"
        )
        args.output_file = f"{tmp_path}/test_filtered.m8"
        args.alignment_fields = ""
        args.filter_field = "alntmscore"
        args.min_val_filter_field = 0.4
        args.max_val_filter_field = 1

        aln_filter_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # First line should be header
        assert lines[0].startswith("query\t")

        # Check that all alignments pass the filter
        header = lines[0].strip().split("\t")
        alntmscore_idx = header.index("alntmscore")

        for line in lines[1:]:
            parts = line.strip().split("\t")
            alntmscore = float(parts[alntmscore_idx])
            assert alntmscore >= args.min_val_filter_field
            assert alntmscore <= args.max_val_filter_field

    def test_filtering_with_explicit_fields(self, tmp_path):
        """Test filtering when alignment_fields are explicitly provided."""
        # Create a test input file without header
        input_file = tmp_path / "input_no_header.m8"
        input_file.write_text(
            "queryA\ttargetA\t0.8\n"
            "queryA\ttargetB\t0.3\n"
            "queryB\ttargetC\t0.5\n"
            "queryB\ttargetD\t0.9\n"
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = "query,target,score"
        args.filter_field = "score"
        args.min_val_filter_field = 0.4
        args.max_val_filter_field = 1.0

        aln_filter_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Should have header + 3 passing alignments (0.8, 0.5, 0.9)
        assert len(lines) == 4
        assert lines[0].strip() == "query\ttarget\tscore"

        # Verify correct alignments kept
        targets = [line.strip().split("\t")[1] for line in lines[1:]]
        assert "targetA" in targets  # 0.8 passes
        assert "targetB" not in targets  # 0.3 fails
        assert "targetC" in targets  # 0.5 passes
        assert "targetD" in targets  # 0.9 passes

    def test_scientific_notation_in_filter_field(self, tmp_path):
        """Test that scientific notation values are properly parsed and filtered."""
        input_file = tmp_path / "input_sci.m8"
        input_file.write_text(
            "query\ttarget\tscore\n"
            "q1\tt1\t4.5E-01\n"
            "q2\tt2\t3.0E-01\n"
            "q3\tt3\t6.0E-01\n"
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.filter_field = "score"
        args.min_val_filter_field = 0.4
        args.max_val_filter_field = 1.0

        aln_filter_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Should have header + 2 passing alignments (0.45 and 0.6)
        assert len(lines) == 3
        # Verify output values are converted to regular float format
        assert "0.45" in lines[1]
        assert "0.6" in lines[2]

    def test_max_less_than_min_raises_error(self, tmp_path):
        """Test that max_val < min_val raises an error."""

        class Args:
            pass

        args = Args()
        args.alignment_file = (
            "tests/test_data/foldseek_related/top_query_per_cluster_tax.m8"
        )
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.filter_field = "alntmscore"
        args.min_val_filter_field = 0.8
        args.max_val_filter_field = 0.4

        with pytest.raises(ValueError, match="max_val can't be less than min_val"):
            aln_filter_main(args)

    def test_missing_filter_field_raises_error(self, tmp_path):
        """Test that a missing filter field raises an error."""
        input_file = tmp_path / "input.m8"
        input_file.write_text("query\ttarget\tscore\nq1\tt1\t0.5\n")

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.filter_field = "nonexistent_field"
        args.min_val_filter_field = 0.0
        args.max_val_filter_field = 1.0

        with pytest.raises(ValueError, match="Cannot find filter field"):
            aln_filter_main(args)

    def test_no_header_no_fields_raises_error(self, tmp_path):
        """Test that missing header without explicit fields raises an error."""
        input_file = tmp_path / "input_no_header.m8"
        input_file.write_text("q1\tt1\t0.5\n")

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.filter_field = "score"
        args.min_val_filter_field = 0.0
        args.max_val_filter_field = 1.0

        with pytest.raises(ValueError, match="alignment_fields has not been provided"):
            aln_filter_main(args)

    def test_all_filtered_out(self, tmp_path):
        """Test behavior when all alignments are filtered out."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tscore\n" "q1\tt1\t0.1\n" "q2\tt2\t0.2\n"
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.filter_field = "score"
        args.min_val_filter_field = 0.5
        args.max_val_filter_field = 1.0

        aln_filter_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Should only have header
        assert len(lines) == 1
        assert lines[0].strip() == "query\ttarget\tscore"

    def test_boundary_values_inclusive(self, tmp_path):
        """Test that min and max boundaries are inclusive."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tscore\n"
            "q1\tt1\t0.4\n"  # exactly at min
            "q2\tt2\t0.6\n"  # in range
            "q3\tt3\t0.8\n"  # exactly at max
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.filter_field = "score"
        args.min_val_filter_field = 0.4
        args.max_val_filter_field = 0.8

        aln_filter_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # All 3 should pass (boundaries are inclusive)
        assert len(lines) == 4  # header + 3 alignments

    def test_skips_duplicate_headers(self, tmp_path):
        """Test that duplicate header lines in the file are skipped."""
        input_file = tmp_path / "input_dup_header.m8"
        input_file.write_text(
            "query\ttarget\tscore\n"
            "q1\tt1\t0.5\n"
            "query\ttarget\tscore\n"  # duplicate header
            "q2\tt2\t0.6\n"
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.alignment_fields = ""
        args.filter_field = "score"
        args.min_val_filter_field = 0.0
        args.max_val_filter_field = 1.0

        aln_filter_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Should have 1 header + 2 alignments
        assert len(lines) == 3
        queries = [line.strip().split("\t")[0] for line in lines]
        assert queries.count("query") == 1  # only one header in output
