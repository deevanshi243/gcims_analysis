# Local maxima peak detection helper for 1D signals.
# Note: GC-IMS data are typically 2D; treat this as a minimal example or extend to 2D.
from typing import Dict, Any
import numpy as np
from scipy.signal import find_peaks

def local_maxima_peaks(signal: np.ndarray, height=None, distance=None, prominence=None, width=None) -> Dict[str, Any]:
    peaks, props = find_peaks(signal, height=height, distance=distance, prominence=prominence, width=width)
    props_serializable = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in props.items()}
    return {"peaks": peaks.tolist(), "properties": props_serializable}
