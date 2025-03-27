# __init__.py
from .solve import (
    solve
)
from .tmm import (
    build_transfer_matrix,
    find_eigenvalues,
)

__all__ = [
    "find_eigenvalues",
    "solve",
    "build_transfer_matrix"
]