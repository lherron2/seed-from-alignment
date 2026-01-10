"""Exact enumeration tests for MCMC sampler correctness.

These tests verify that the MCMC sampler produces the correct
Boltzmann distribution by comparing against exact enumeration
for small systems.

Acceptance criteria: max frequency error < 5%
"""

import pytest

# Tests will be implemented in Phase 2


class TestExactDistribution:
    """Tests comparing MCMC samples to exact Boltzmann distribution."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_tiny_system_toggle_only(self) -> None:
        """
        Test on tiny system (L=8, ~10 candidates) with toggle moves only.

        Procedure:
        1. Enumerate all valid matchings
        2. Compute exact Boltzmann probabilities
        3. Run MCMC for 10,000 samples
        4. Compare empirical frequencies to exact probabilities
        5. Assert max frequency error < 5%
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_tiny_system_with_segments(self) -> None:
        """
        Test on tiny system with segment moves enabled.

        Same procedure as toggle-only, but with segment birth/death.
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 4")
    def test_tiny_system_full_moves(self) -> None:
        """
        Test on tiny system with all move types enabled.

        Same procedure, but includes loop closure moves.
        """
        pass


class TestDetailedBalance:
    """Tests for detailed balance verification."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_toggle_detailed_balance(self) -> None:
        """
        Verify π(M)·A(M→M') = π(M')·A(M'→M) for toggle moves.

        For many random (M, M') pairs:
        1. Compute acceptance probabilities both directions
        2. Compute Boltzmann weights
        3. Verify equation within tolerance (1e-8)
        """
        pass

    @pytest.mark.skip(reason="Not implemented yet - Phase 3")
    def test_segment_detailed_balance(self) -> None:
        """
        Verify detailed balance for segment moves.

        Must account for Hastings ratios.
        """
        pass


class TestErgodicity:
    """Tests for MCMC ergodicity."""

    @pytest.mark.skip(reason="Not implemented yet - Phase 2")
    def test_all_states_reachable(self) -> None:
        """
        For tiny system, verify all valid matchings are sampled.
        """
        pass
