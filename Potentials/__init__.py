# Explicitly expose functions/classes for cleaner imports
from .pot import *
from .discret import *

# Optional: Shorten imports for users
__all__ = ["pot","pot_period","discretiser_period", "discretiser"]