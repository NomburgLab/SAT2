from sat.scripts.struc_get_contact_probability import (
    sum_contact_probabilities,
    highest_contact_probability,
)

import pytest
import numpy as np
import scipy.special


class TestSumContactProbabilities:
    """Test suite for the sum_contact_probabilities function."""

    @pytest.fixture
    def example_distogram(self):
        """Fixture to provide a sample distogram for testing."""
        return {
            "logits": np.random.rand(5, 5, 64),  # Example logits for a small test case
            "bin_edges": np.linspace(2.3, 22.8, 63),  # Example bin edges
        }

    def test_basic_functionality(self, example_distogram):
        """Test basic functionality with a valid length cutoff."""
        length_cutoff = 8.0
        result = sum_contact_probabilities(example_distogram, length_cutoff)

        assert result.shape == (5, 5)
        assert np.all(result >= 0), "Contact probabilities should be non-negative."

    def test_length_cutoff(self, example_distogram):
        """Ensure probabilities are summed correctly for bins ≤ length_cutoff."""
        length_cutoff = 10.0
        logits = example_distogram["logits"]
        probabilities = scipy.special.softmax(logits, axis=-1)

        # Manually compute expected results for bins ≤ length_cutoff
        valid_bins = np.where(example_distogram["bin_edges"] <= length_cutoff)[0]
        expected = np.sum(probabilities[:, :, valid_bins], axis=-1)

        result = sum_contact_probabilities(example_distogram, length_cutoff)
        np.testing.assert_array_almost_equal(result, expected, decimal=5)

    def test_empty_distogram(self):
        """Test behavior when distogram is missing required keys."""
        incomplete_distogram = {
            "logits": np.random.rand(5, 5, 64)
        }  # Missing "bin_edges"
        with pytest.raises(
            ValueError, match="Distogram must contain 'logits' and 'bin_edges'."
        ):
            sum_contact_probabilities(incomplete_distogram, 8.0)

    def test_no_bins_within_cutoff(self, example_distogram):
        """Ensure the function raises an error if no bins satisfy length_cutoff."""
        length_cutoff = 1.0  # Smaller than the smallest bin edge
        with pytest.raises(
            ValueError, match="No bins found for the given length cutoff"
        ):
            sum_contact_probabilities(example_distogram, length_cutoff)

    def test_output_shape(self, example_distogram):
        """Check that the output matrix has the correct shape."""
        length_cutoff = 12.0
        logits_shape = example_distogram["logits"].shape
        result = sum_contact_probabilities(example_distogram, length_cutoff)

        assert result.shape == logits_shape[:2], "Output shape should match (L, L)."

    @pytest.fixture
    def small_distogram(self):
        """Fixture providing a small, readable distogram."""
        return {
            "logits": np.array(
                [[[1, 2, 0, 0], [0, 1, 0, 0]], [[2, 1, 0, 0], [1, 0, 0, 0]]]
            ),  # 2x2 logits with 4 bins
            "bin_edges": np.array([2.0, 5.0, 10.0, 15.0]),  # 4 bins
        }

    def test_small_distogram_with_cutoff_5(self, small_distogram):
        """Test with a small distogram and cutoff of 5.0 Å."""
        length_cutoff = 5.0
        result = sum_contact_probabilities(small_distogram, length_cutoff)

        # Expected output
        # Softmax probabilities along last axis
        probabilities = scipy.special.softmax(small_distogram["logits"], axis=-1)
        # Include bins corresponding to distances ≤ 5.0
        expected = np.sum(probabilities[:, :, :2], axis=-1)  # Use first two bins

        np.testing.assert_array_almost_equal(result, expected, decimal=5)

    def test_small_distogram_with_cutoff_10(self, small_distogram):
        """Test with a small distogram and cutoff of 10.0 Å."""
        length_cutoff = 10.0
        result = sum_contact_probabilities(small_distogram, length_cutoff)

        # Expected output
        # Softmax probabilities along last axis
        probabilities = scipy.special.softmax(small_distogram["logits"], axis=-1)
        # Include bins corresponding to distances ≤ 10.0
        expected = np.sum(probabilities[:, :, :3], axis=-1)  # Use first three bins

        np.testing.assert_array_almost_equal(result, expected, decimal=5)


class TestHighestContactProbability:
    """Test suite for the highest_contact_probability function using small arrays."""

    def test_highest_contact_probability_basic(self):
        """Test with a simple 2x2 contact probability matrix."""
        # Small readable summed probabilities matrix
        summed_probabilities = np.array(
            [
                [0.1, 0.2, 0.3, 0.4],
                [0.5, 0.6, 0.7, 0.8],
                [0.9, 1.0, 1.1, 1.2],
                [1.3, 1.4, 1.5, 1.6],
            ]
        )
        len1 = 2  # Protein 1 length
        # Submatrix: summed_probabilities[:2, 2:4] = [[0.3, 0.4], [0.7, 0.8]]
        expected = 0.8  # Highest value in submatrix

        result = highest_contact_probability(summed_probabilities, len1)
        assert result == expected, f"Expected {expected}, but got {result}."

    def test_highest_contact_probability_single_row(self):
        """Test when Protein 1 consists of a single residue."""
        # Small readable summed probabilities matrix
        summed_probabilities = np.array(
            [
                [0.2, 0.3, 0.4, 0.5],
                [0.6, 0.7, 0.8, 0.9],
                [1.0, 1.1, 1.2, 1.3],
                [1.4, 1.5, 1.6, 1.7],
            ]
        )
        len1 = 1  # Protein 1 length
        # Submatrix: summed_probabilities[:1, 1:] = [[0.3, 0.4, 0.5]]
        expected = 0.5  # Highest value in submatrix

        result = highest_contact_probability(summed_probabilities, len1)
        assert result == expected, f"Expected {expected}, but got {result}."

    def test_highest_contact_probability_equal_values(self):
        """Test when all values in the submatrix are the same."""
        # Small readable summed probabilities matrix
        summed_probabilities = np.array(
            [
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1],
            ]
        )
        len1 = 2  # Protein 1 length
        # Submatrix: summed_probabilities[:2, 2:4] = [[0.1, 0.1], [0.1, 0.1]]
        expected = 0.1  # All values are equal

        result = highest_contact_probability(summed_probabilities, len1)
        assert result == expected, f"Expected {expected}, but got {result}."

    def test_highest_contact_probability_invalid_lengths(self):
        """Test when Protein 1 length is invalid."""
        # Small readable summed probabilities matrix
        summed_probabilities = np.array(
            [[0.2, 0.3, 0.4], [0.5, 0.6, 0.7], [0.8, 0.9, 1.0]]
        )
        len1 = 4  # Invalid length (exceeds matrix dimensions)

        with pytest.raises(ValueError, match="Invalid Protein 1 length"):
            highest_contact_probability(summed_probabilities, len1)
