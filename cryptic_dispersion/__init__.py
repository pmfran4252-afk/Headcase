"""NMR-inspired detection of cryptic pockets from molecular-dynamics ensembles.

The idea in one line: relaxation dispersion finds cryptic pockets because it
measures exchange with a sparsely populated excited state, so reproduce that
measurement on an MD ensemble -- but with an observable chosen to be sensitive
to pocket opening rather than one that responds to it by accident, and with the
structure of the excited state in hand rather than invisible.
"""

from .exchange import (
    ExchangeFit,
    detectability_limit,
    dispersion_amplitude,
    exchange_contrast,
    fit_two_state,
    is_detectable,
    nmr_window,
    r2eff_carver_richards,
    r2eff_from_fit,
    r2eff_luz_meiboom,
)
from .hdx import AmideEnvironment, amide_environment, breathing_anomaly, ln_protection_factor
from .observables import Projection, contact_number, project_all, robust_z, slow_projection
from .score import SiteScore, analyse, rank_sites
from .tails import TailFit, fit_tail, latent_openness, threshold_stability

__version__ = "0.1.0"

__all__ = [
    "ExchangeFit", "fit_two_state", "exchange_contrast", "dispersion_amplitude",
    "detectability_limit", "is_detectable", "nmr_window",
    "r2eff_carver_richards", "r2eff_luz_meiboom", "r2eff_from_fit",
    "TailFit", "fit_tail", "latent_openness", "threshold_stability",
    "AmideEnvironment", "amide_environment", "breathing_anomaly", "ln_protection_factor",
    "Projection", "slow_projection", "project_all", "contact_number", "robust_z",
    "SiteScore", "analyse", "rank_sites",
    "__version__",
]
