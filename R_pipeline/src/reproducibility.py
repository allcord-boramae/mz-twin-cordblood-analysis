import os
import random
import numpy as np
from config import GLOBAL_SEED

def seed_everything(seed: int = GLOBAL_SEED) -> int:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    return seed
