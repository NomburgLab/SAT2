import pytest
from sat.scripts.aln_cigar_to_cov import (
    aln_cigar_to_cov_main,
    parse_cigar_matches,
    parse_alignment_file,
    add_cigar_coverage_columns,
)


class TestParseCigarMatches:
    """Tests for the parse_cigar_matches helper function."""

    def test_simple_match(self):
        """Test parsing a simple all-match CIGAR string."""
        assert parse_cigar_matches("118M") == 118
        assert parse_cigar_matches("50M") == 50

    def test_multiple_matches(self):
        """Test parsing CIGAR with multiple M operations."""
        # 6M + 19M + 1M + 1M + 19M + 8M + 30M = 84
        assert parse_cigar_matches("6M9I19M3I1M1I1M1I19M1D8M1D30M") == 84

    def test_complex_cigar(self):
        """Test a complex CIGAR string from real foldseek output."""
        # This is from the test data file, row 2
        # 1M26I3M2I1M1D10I5M1D2I6M1D7I10M1I10M1D1M1D5M1I11M84D3I13M
        # Counts: 1+3+1+5+6+10+10+1+5+11+13 = 66
        cigar = "1M26I3M2I1M1D10I5M1D2I6M1D7I10M1I10M1D1M1D5M1I11M84D3I13M"
        assert parse_cigar_matches(cigar) == 66

    def test_no_matches(self):
        """Test CIGAR with no M operations."""
        assert parse_cigar_matches("10I5D") == 0

    def test_empty_string(self):
        """Test empty CIGAR string."""
        assert parse_cigar_matches("") == 0


class TestParseAlignmentFile:
    """Tests for file parsing."""

    def test_parse_with_auto_header(self, tmp_path):
        """Test parsing file with automatic header detection."""
        input_file = tmp_path / "test.m8"
        input_file.write_text(
            "query\ttarget\tqlen\ttlen\tcigar\n"
            "q1\tt1\t100\t100\t100M\n"
            "q2\tt2\t50\t50\t50M\n"
        )

        columns, data_rows = parse_alignment_file(str(input_file))

        assert columns == ["query", "target", "qlen", "tlen", "cigar"]
        assert len(data_rows) == 2
        assert data_rows[0][0] == "q1"
        assert data_rows[1][0] == "q2"

    def test_parse_with_provided_colnames(self, tmp_path):
        """Test parsing file with explicitly provided column names."""
        input_file = tmp_path / "test_no_header.m8"
        input_file.write_text(
            "q1\tt1\t100\t100\t100M\n"
            "q2\tt2\t50\t50\t50M\n"
        )

        colnames = "query,target,qlen,tlen,cigar"
        columns, data_rows = parse_alignment_file(
            str(input_file), colnames
        )

        assert columns == ["query", "target", "qlen", "tlen", "cigar"]
        assert len(data_rows) == 2

    def test_parse_with_provided_colnames_and_header(self, tmp_path):
        """Test that header is skipped when colnames provided and header exists."""
        input_file = tmp_path / "test.m8"
        input_file.write_text(
            "query\ttarget\tqlen\ttlen\tcigar\n"
            "q1\tt1\t100\t100\t100M\n"
        )

        colnames = "query,target,qlen,tlen,cigar"
        columns, data_rows = parse_alignment_file(
            str(input_file), colnames
        )

        assert len(data_rows) == 1  # Header skipped

    def test_missing_required_column_raises_error(self, tmp_path):
        """Test that missing required column raises error."""
        input_file = tmp_path / "test.m8"
        input_file.write_text(
            "query\ttarget\tqlen\tcigar\n"  # missing tlen
            "q1\tt1\t100\t100M\n"
        )

        with pytest.raises(ValueError, match="Required column 'tlen' not found"):
            parse_alignment_file(str(input_file))

    def test_no_header_no_colnames_raises_error(self, tmp_path):
        """Test that missing header without colnames raises error."""
        input_file = tmp_path / "test.m8"
        input_file.write_text("q1\tt1\t100\t100\t100M\n")

        with pytest.raises(ValueError, match="colnames not provided"):
            parse_alignment_file(str(input_file))


class TestAddCigarCoverageColumns:
    """Tests for adding coverage columns."""

    def test_columns_inserted_after_cigar(self):
        """Test that cigar_qcov and cigar_tcov are inserted after cigar column."""
        columns = ["query", "target", "qlen", "tlen", "cigar", "evalue"]
        data_rows = [
            ["q1", "t1", "100", "100", "80M", "1e-10"],
        ]

        new_columns, new_rows = add_cigar_coverage_columns(columns, data_rows)

        assert new_columns == [
            "query", "target", "qlen", "tlen", "cigar",
            "cigar_qcov", "cigar_tcov", "evalue"
        ]
        assert len(new_rows) == 1
        assert new_rows[0][5] == "0.800"  # cigar_qcov
        assert new_rows[0][6] == "0.800"  # cigar_tcov
        assert new_rows[0][7] == "1e-10"  # evalue preserved

    def test_coverage_calculation(self):
        """Test correct coverage calculation."""
        columns = ["query", "target", "qlen", "tlen", "cigar"]
        data_rows = [
            ["q1", "t1", "100", "200", "80M"],  # qcov=0.8, tcov=0.4
            ["q2", "t2", "50", "100", "25M"],   # qcov=0.5, tcov=0.25
        ]

        new_columns, new_rows = add_cigar_coverage_columns(columns, data_rows)

        assert new_rows[0][5] == "0.800"  # cigar_qcov
        assert new_rows[0][6] == "0.400"  # cigar_tcov
        assert new_rows[1][5] == "0.500"  # cigar_qcov
        assert new_rows[1][6] == "0.250"  # cigar_tcov

    def test_zero_length_handled(self):
        """Test that zero length sequences get 0.0 coverage."""
        columns = ["query", "target", "qlen", "tlen", "cigar"]
        data_rows = [
            ["q1", "t1", "0", "100", "80M"],   # zero qlen
            ["q2", "t2", "100", "0", "80M"],   # zero tlen
        ]

        new_columns, new_rows = add_cigar_coverage_columns(columns, data_rows)

        assert new_rows[0][5] == "0.000"  # cigar_qcov (zero qlen)
        assert new_rows[0][6] == "0.800"  # cigar_tcov
        assert new_rows[1][5] == "0.800"  # cigar_qcov
        assert new_rows[1][6] == "0.000"  # cigar_tcov (zero tlen)

    def test_cigar_at_end_of_columns(self):
        """Test when cigar is the last column."""
        columns = ["query", "target", "qlen", "tlen", "cigar"]
        data_rows = [
            ["q1", "t1", "100", "100", "100M"],
        ]

        new_columns, new_rows = add_cigar_coverage_columns(columns, data_rows)

        assert new_columns == [
            "query", "target", "qlen", "tlen", "cigar",
            "cigar_qcov", "cigar_tcov"
        ]
        assert new_rows[0] == ["q1", "t1", "100", "100", "100M", "1.000", "1.000"]


class TestAlnCigarToCovMain:
    """Integration tests for the main function."""

    def test_basic_with_real_data(self, tmp_path):
        """Test with real test data."""

        class Args:
            pass

        args = Args()
        args.alignment_file = (
            "tests/test_data/foldseek_related/cigar_filter_test.m8"
        )
        args.output_file = f"{tmp_path}/output.m8"
        args.colnames = ""

        aln_cigar_to_cov_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Should have header + all data rows
        assert len(lines) >= 1
        assert lines[0].startswith("query\t")

        # Check that cigar_qcov and cigar_tcov columns are added after cigar
        header = lines[0].strip().split("\t")
        assert "cigar_qcov" in header
        assert "cigar_tcov" in header

        cigar_idx = header.index("cigar")
        assert header[cigar_idx + 1] == "cigar_qcov"
        assert header[cigar_idx + 2] == "cigar_tcov"

    def test_columns_positioned_correctly(self, tmp_path):
        """Test that new columns are placed after cigar."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tqlen\ttlen\tcigar\tevalue\tbits\n"
            "q1\tt1\t100\t100\t100M\t1e-10\t200\n"
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.colnames = ""

        aln_cigar_to_cov_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        header = lines[0].strip().split("\t")
        assert header == [
            "query", "target", "qlen", "tlen", "cigar",
            "cigar_qcov", "cigar_tcov", "evalue", "bits"
        ]

        # Check data row
        data = lines[1].strip().split("\t")
        assert data[5] == "1.000"  # cigar_qcov
        assert data[6] == "1.000"  # cigar_tcov
        assert data[7] == "1e-10"  # evalue preserved
        assert data[8] == "200"    # bits preserved

    def test_all_rows_processed(self, tmp_path):
        """Test that all rows are output (no filtering)."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tqlen\ttlen\tcigar\n"
            "q1\tt1\t100\t100\t90M\n"
            "q2\tt2\t100\t100\t50M\n"
            "q3\tt3\t100\t100\t10M\n"
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.colnames = ""

        aln_cigar_to_cov_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # All 3 data rows should be present
        assert len(lines) == 4  # header + 3 data rows

    def test_explicit_colnames(self, tmp_path):
        """Test with explicitly provided column names."""
        input_file = tmp_path / "input_no_header.m8"
        input_file.write_text(
            "q1\tt1\t100\t100\t100M\n"
            "q2\tt2\t100\t100\t40M\n"
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.colnames = "query,target,qlen,tlen,cigar"

        aln_cigar_to_cov_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        # Header should include new columns
        header = lines[0].strip().split("\t")
        assert header == [
            "query", "target", "qlen", "tlen", "cigar",
            "cigar_qcov", "cigar_tcov"
        ]
        assert len(lines) == 3  # header + 2 data rows

    def test_coverage_values_correct(self, tmp_path):
        """Test that coverage values are calculated correctly."""
        input_file = tmp_path / "input.m8"
        input_file.write_text(
            "query\ttarget\tqlen\ttlen\tcigar\n"
            "q1\tt1\t100\t200\t80M\n"  # qcov=0.8, tcov=0.4
        )

        class Args:
            pass

        args = Args()
        args.alignment_file = str(input_file)
        args.output_file = f"{tmp_path}/output.m8"
        args.colnames = ""

        aln_cigar_to_cov_main(args)

        with open(args.output_file) as f:
            lines = f.readlines()

        data = lines[1].strip().split("\t")
        assert data[5] == "0.800"  # cigar_qcov
        assert data[6] == "0.400"  # cigar_tcov
