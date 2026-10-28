"""Structure generators for the discovery loop (DISCOVERY.md 2.2, 4).

``rtsc.generate.mattergen``: Microsoft MatterGen through its CLI in a separate environment.
``rtsc.generate.finetune_export``: ledger labels -> MatterGen fine-tuning CSV dataset.
Every generated structure carries the HYPOTHETICAL_STRUCTURE flag.  Nothing here imports
MatterGen or torch.
"""
from .mattergen import (CHECKPOINT_PROPERTIES, MatterGenConfig, MatterGenGenerator, MatterGenRunError,
                        UnsupportedConditioningError)

__all__ = ["MatterGenGenerator", "MatterGenConfig", "UnsupportedConditioningError", "MatterGenRunError",
           "CHECKPOINT_PROPERTIES"]
