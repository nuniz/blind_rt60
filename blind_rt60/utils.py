import numpy as np


def calculate_decay_time(decay_db: float, tau: float) -> float:
    """
    Calculates the decay time of a signal based on its decay in decibels (dB) and the time constant (tau).
    decay_time = -decay_db / (20 * log10(e)) * tau

    Parameters:
        - decay_db (float): The decay of the signal in decibels (dB). Positive values indicate
                             attenuation, while negative values indicate amplification.
        - tau (float): The time constant of the system, which represents the time it takes for
                       the signal to decay to 1/e (approximately 36.8%) of its initial value.
        - e (float, optional): The mathematical constant e (approximately 2.71828). Defaults
                                to `np.e` for efficiency (already imported with `numpy`).

    Returns:
        float: The calculated decay time in the same units as `tau` (typically seconds, milliseconds, etc.).

    Raises:
        ValueError: If `decay_db` is not a finite number (i.e., NaN or Inf).
    """

    if not np.isfinite(decay_db):
        raise ValueError("decay_db must be a finite number (not NaN or Inf).")

    decay_time = -decay_db / (20 * np.log10(np.e**-1)) * tau
    return decay_time
