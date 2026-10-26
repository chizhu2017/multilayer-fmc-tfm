from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class Config:
    nx: int = 600
    h1: int = 350
    h2: int = 200
    h3: int = 50
    v1: float = 4200.0
    v2: float = 2800.0
    v3: float = 4600.0
    x_loc: int = 50
    gap: int = 50
    cut_point: float = 20.0
    t0_fine: float = 0.0
    color_ratio: float = 0.7
    show_interface: bool = False

    def __post_init__(self):
        for key in ("nx", "h1", "h2", "h3", "x_loc", "gap"):
            value = getattr(self, key)
            if type(value) is not int or value <= 0:
                raise ValueError(f"{key} must be a positive integer")
        if self.h1 < 2:
            raise ValueError("h1 must be >= 2 mm")
        for key in ("v1", "v2", "v3", "color_ratio"):
            value = getattr(self, key)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{key} must be finite and positive")
        for key in ("cut_point", "t0_fine"):
            if not math.isfinite(getattr(self, key)):
                raise ValueError(f"{key} must be finite")
        if type(self.show_interface) is not bool:
            raise ValueError("show_interface must be boolean")

    @property
    def nz(self):
        return self.h1 + self.h2 + self.h3

    @classmethod
    def load(cls, path):
        return cls(**json.loads(Path(path).read_text(encoding="utf-8")))

    def to_dict(self):
        return asdict(self)
