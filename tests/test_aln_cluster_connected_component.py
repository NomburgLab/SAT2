import pytest
import tempfile
import os

from sat.scripts.aln_cluster_connected_component import (
    find_connected_components,
    parse_alignment_file,
)


class TestFindConnectedComponents:
    """Tests for the find_connected_components function."""

    def test_self_alignments_only_separate_clusters(self):
        """
        Proteins that only align to themselves should be in separate clusters.
        """
        # Each protein only has a self-alignment
        pairs = [
            ("A", "A"),
            ("B", "B"),
            ("C", "C"),
        ]

        clusters = find_connected_components(pairs)

        # Should have 3 separate single-member clusters
        assert len(clusters) == 3

        # Convert to sorted list of sorted members for comparison
        cluster_members = sorted([sorted(list(c)) for c in clusters])
        assert cluster_members == [["A"], ["B"], ["C"]]

    def test_self_alignments_mixed_with_real_alignments(self):
        """
        Self-alignments mixed with real alignments.
        A and B only self-align, C and D align to each other.
        """
        pairs = [
            ("A", "A"),
            ("B", "B"),
            ("C", "C"),
            ("D", "D"),
            ("C", "D"),  # C and D are connected
        ]

        clusters = find_connected_components(pairs)

        # Should have 3 clusters: {A}, {B}, {C, D}
        assert len(clusters) == 3

        cluster_members = sorted([sorted(list(c)) for c in clusters])
        assert cluster_members == [["A"], ["B"], ["C", "D"]]

    def test_transitive_connections(self):
        """
        Test that transitive connections work correctly.
        A->B, B->C should put A, B, C in the same cluster.
        """
        pairs = [
            ("A", "B"),
            ("B", "C"),
        ]

        clusters = find_connected_components(pairs)

        # Should have 1 cluster containing A, B, C
        assert len(clusters) == 1
        assert clusters[0] == {"A", "B", "C"}

    def test_two_separate_clusters(self):
        """
        Two separate clusters that don't connect.
        """
        pairs = [
            ("A", "B"),
            ("C", "D"),
        ]

        clusters = find_connected_components(pairs)

        # Should have 2 clusters
        assert len(clusters) == 2

        cluster_members = sorted([sorted(list(c)) for c in clusters])
        assert cluster_members == [["A", "B"], ["C", "D"]]

    def test_single_connector_merges_clusters(self):
        """
        A single alignment can merge two otherwise separate groups.
        This demonstrates how one bad/promiscuous hit can merge everything.
        """
        pairs = [
            # Group 1
            ("A", "B"),
            ("B", "C"),
            # Group 2
            ("X", "Y"),
            ("Y", "Z"),
            # Connector - this single alignment merges everything!
            ("C", "X"),
        ]

        clusters = find_connected_components(pairs)

        # Should have 1 cluster containing all members
        assert len(clusters) == 1
        assert clusters[0] == {"A", "B", "C", "X", "Y", "Z"}

    def test_empty_pairs(self):
        """
        Empty pairs list should return empty clusters.
        """
        pairs = []
        clusters = find_connected_components(pairs)
        assert len(clusters) == 0

    def test_all_inputs_adds_missing_members(self):
        """
        Test that all_inputs adds members not in pairs as single-member clusters.
        """
        pairs = [
            ("A", "B"),
        ]
        all_inputs = {"A", "B", "C", "D"}

        clusters = find_connected_components(pairs, all_inputs)

        # Should have 3 clusters: {A, B}, {C}, {D}
        assert len(clusters) == 3

        cluster_members = sorted([sorted(list(c)) for c in clusters])
        assert cluster_members == [["A", "B"], ["C"], ["D"]]

    def test_all_inputs_with_self_alignments(self):
        """
        Test all_inputs with self-alignments.
        Members that self-align should NOT be duplicated.
        """
        pairs = [
            ("A", "A"),
            ("B", "B"),
        ]
        all_inputs = {"A", "B", "C"}

        clusters = find_connected_components(pairs, all_inputs)

        # Should have 3 clusters: {A}, {B}, {C}
        assert len(clusters) == 3

        cluster_members = sorted([sorted(list(c)) for c in clusters])
        assert cluster_members == [["A"], ["B"], ["C"]]

    def test_duplicate_pairs_dont_affect_result(self):
        """
        Duplicate pairs shouldn't affect the clustering.
        """
        pairs = [
            ("A", "B"),
            ("A", "B"),  # duplicate
            ("B", "A"),  # reverse
            ("C", "C"),
        ]

        clusters = find_connected_components(pairs)

        # Should have 2 clusters: {A, B}, {C}
        assert len(clusters) == 2

        cluster_members = sorted([sorted(list(c)) for c in clusters])
        assert cluster_members == [["A", "B"], ["C"]]

    def test_large_chain_stays_connected(self):
        """
        A long chain should stay as one cluster.
        """
        pairs = [
            ("A", "B"),
            ("B", "C"),
            ("C", "D"),
            ("D", "E"),
            ("E", "F"),
        ]

        clusters = find_connected_components(pairs)

        assert len(clusters) == 1
        assert clusters[0] == {"A", "B", "C", "D", "E", "F"}

    def test_star_topology(self):
        """
        A star topology where one node connects to many others.
        """
        pairs = [
            ("hub", "A"),
            ("hub", "B"),
            ("hub", "C"),
            ("hub", "D"),
        ]

        clusters = find_connected_components(pairs)

        assert len(clusters) == 1
        assert clusters[0] == {"hub", "A", "B", "C", "D"}


class TestParseAlignmentFile:
    """Tests for parsing alignment files."""

    def test_parse_with_header(self):
        """Test parsing a file with a header."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".m8", delete=False) as f:
            f.write("query\ttarget\tscore\n")
            f.write("A\tA\t1.0\n")
            f.write("B\tB\t1.0\n")
            f.write("C\tD\t0.8\n")
            f.name

        try:
            pairs = parse_alignment_file(f.name)
            assert pairs == [("A", "A"), ("B", "B"), ("C", "D")]
        finally:
            os.unlink(f.name)

    def test_parse_with_provided_colnames(self):
        """Test parsing a file with provided column names."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".m8", delete=False) as f:
            f.write("A\tA\t1.0\n")
            f.write("B\tB\t1.0\n")
            f.write("C\tD\t0.8\n")
            f.name

        try:
            pairs = parse_alignment_file(f.name, colnames="query,target,score")
            assert pairs == [("A", "A"), ("B", "B"), ("C", "D")]
        finally:
            os.unlink(f.name)


class TestIntegration:
    """Integration tests simulating real-world scenarios."""

    def test_scenario_self_alignments_become_giant_cluster(self):
        """
        Reproduce the reported bug: proteins with only self-alignments
        ending up in one giant cluster.

        This test verifies that self-alignments alone do NOT create
        connections between different proteins.
        """
        # Simulate an alignment file where each protein only aligns to itself
        pairs = [
            ("protein1", "protein1"),
            ("protein2", "protein2"),
            ("protein3", "protein3"),
            ("protein4", "protein4"),
            ("protein5", "protein5"),
        ]

        clusters = find_connected_components(pairs)

        # CRITICAL: Each protein should be in its own cluster
        assert len(clusters) == 5, (
            f"Expected 5 separate clusters for 5 proteins with only self-alignments, "
            f"but got {len(clusters)} clusters"
        )

        # Verify each cluster has exactly one member
        for cluster in clusters:
            assert len(cluster) == 1, (
                f"Expected each cluster to have 1 member, "
                f"but found cluster with {len(cluster)} members: {cluster}"
            )

    def test_scenario_one_promiscuous_hit_merges_all(self):
        """
        Demonstrate how one promiscuous/spurious alignment can merge
        otherwise separate clusters.

        This is the expected behavior - it's not a bug, but it shows
        why strict filtering is needed before clustering.
        """
        # 5 proteins, each only self-aligning
        pairs = [
            ("protein1", "protein1"),
            ("protein2", "protein2"),
            ("protein3", "protein3"),
            ("protein4", "protein4"),
            ("protein5", "protein5"),
        ]

        clusters_before = find_connected_components(pairs)
        assert len(clusters_before) == 5, "Should have 5 separate clusters"

        # Now add ONE spurious alignment that connects protein1 to protein2
        pairs_with_spurious = pairs + [("protein1", "protein2")]

        clusters_after = find_connected_components(pairs_with_spurious)

        # Now protein1 and protein2 should be merged
        assert len(clusters_after) == 4, (
            "After adding one connection, should have 4 clusters"
        )

    def test_debug_real_data_pattern(self):
        """
        Debug test to understand what's happening with real data.
        If you're seeing all proteins in one cluster, check if there
        are any cross-alignments in your data.
        """
        # Simulate what might happen in real data
        pairs = [
            # Self alignments
            ("A", "A"),
            ("B", "B"),
            ("C", "C"),
            # Maybe there's a hidden alignment somewhere?
        ]

        clusters = find_connected_components(pairs)

        print("\n=== DEBUG OUTPUT ===")
        print(f"Input pairs: {pairs}")
        print(f"Number of clusters: {len(clusters)}")
        for i, cluster in enumerate(clusters):
            print(f"  Cluster {i}: {cluster}")
        print("===================\n")

        # This should pass if the algorithm is correct
        assert len(clusters) == 3


