# __init__.py
from .solve import (
    solve,
    solve_discretized
)
from .tmm import (
    build_transfer_matrix,
    find_eigenvalues,
    k_from_energy,
)

__all__ = [
    "find_eigenvalues",
    "solve",
    "build_transfer_matrix",
    "k_from_energy",
    "solve_discretized"]