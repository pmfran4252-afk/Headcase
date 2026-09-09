"""What this trajectory could have seen, separately from what it did see.

A residue that scores low has told you one of two completely different things:
that it was watched and stayed shut, or that nothing was ever going to be
visible at this trajectory length.  Collapsing both into a low score reports a
ranking conditional on being observable, which is not the question anyone asked.

At ATLAS lengths the second case is the common one.  Cryptic opening runs at
roughly 1e3-1e6 s^-1; a 100 ns trajectory needs about 1e9 s^-1 before a
two-state rate can be fitted at all.  So the useful output is not a number per
residue but a per-channel verdict:

``seen``
    the channel produced a usable value.

``no_signal``
    the channel ran, the floor says it could have resolved a reference site,
    and it found nothing.  This is an informative negative.

``blind``
    the channel could not have resolved a reference site at this trajectory
    length.  Its silence carries no information at all.

and a residue-level outcome of ``SCORED``, ``NEGATIVE`` or ``ABSENT``.

The floors are arithmetic, not thresholds chosen to taste.  The expected number
of A->B transitions in a trajectory of length T is ``T p_A p_B k_ex``, so any
method needing ``m`` of them cannot resolve exchange slower than
``m / (T p_A p_B)``.  That single relation gives every floor below.

One consequence is worth stating plainly because it inverts the documented
design.  The dispersion channel needs ``min_transitions`` (8) transitions; the
tail channel needs ``min_clusters`` (25) independent excursions, and
declustering merges everything inside one autocorrelation time, so its
excursion count can never exceed the number of open/close cycles.  The tail
channel therefore needs an exchange rate about 6x *faster* than dispersion --
the opposite of surviving the undersampling that silences dispersion.  What
would restore the intended behaviour is counting partial excursions rather than
declustering them away; :func:`tail_floor` states the cost of not doing so.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Optional, Sequence

import numpy as np

from .tails import autocorr_time, fit_tail

# Per-channel verdicts.
SEEN = "seen"
NO_SIGNAL = "no_signal"
BLIND = "blind"

# Residue-level outcomes.
SCORED = "SCORED"
NEGATIVE = "NEGATIVE"
ABSENT = "ABSENT"

# Exchange rate a genuine cryptic opening runs at, against which every floor is
# judged.  Cryptic-pocket opening is commonly 1e3-1e6 s^-1; 1e6 is the generous
# end, so a channel blind at 1e6 is blind at every rate that matters.
REFERENCE_K_EX = 1.0e6
REFERENCE_POPULATION = 0.05


@dataclass
class ChannelStatus:
    channel: str
    status: str
    floor: float
    floor_units: str
    note: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class ResidueStatus:
    residue: int
    outcome: str
    channels: dict = field(default_factory=dict)

    @property
    def blind_channels(self) -> list[str]:
        return [k for k, v in self.channels.items() if v.status == BLIND]

    def as_dict(self) -> dict:
        return {"residue": self.residue, "outcome": self.outcome,
                "channels": {k: v.as_dict() for k, v in self.channels.items()}}


# --------------------------------------------------------------------------
# Floors.  All three follow from E[transitions] = T p_A p_B k_ex.
# --------------------------------------------------------------------------

def dispersion_floor(total_time_ns: float, p_b: float = REFERENCE_POPULATION,
                     min_transitions: int = 8) -> float:
    """Slowest ``k_ex`` (s^-1) a two-state fit can resolve in this trajectory."""
    if not 0.0 < p_b < 1.0 or total_time_ns <= 0:
        return float("inf")
    return float(min_transitions / (2.0 * total_time_ns * 1e-9 * p_b * (1.0 - p_b)))


def tail_floor(total_time_ns: float, p_b: float = REFERENCE_POPULATION,
               min_clusters: int = 25) -> float:
    """Slowest ``k_ex`` the declustered tail fit can support.

    Declustering with a gap of one autocorrelation time merges every excursion
    inside a single opening episode, so the independent-excursion count is
    bounded above by the number of open/close cycles.  Requiring
    ``min_clusters`` of them is therefore a rate floor, and a stricter one than
    the dispersion channel's.
    """
    if not 0.0 < p_b < 1.0 or total_time_ns <= 0:
        return float("inf")
    return float(min_clusters / (total_time_ns * 1e-9 * p_b * (1.0 - p_b)))


def breathing_floor(n_frames: int, tau_frames: float = 1.0,
                    min_open_frames: int = 10) -> float:
    """Rarest state (as a population) whose exchange the ensemble average reflects.

    A population floor rather than a rate floor: the protection factor is a time
    average, so it needs occupied frames, not completed round trips.  That is
    why this channel survives where the other two go blind.
    """
    if n_frames <= 0:
        return 1.0
    return float(min(1.0, min_open_frames * max(tau_frames, 1.0) / n_frames))


# --------------------------------------------------------------------------
# Per-channel status.
# --------------------------------------------------------------------------

def dispersion_status(fit, total_time_ns: float, *,
                      reference_k_ex: float = REFERENCE_K_EX,
                      reference_population: float = REFERENCE_POPULATION,
                      min_transitions: int = 8) -> ChannelStatus:
    floor = dispersion_floor(total_time_ns, reference_population, min_transitions)
    if floor > reference_k_ex:
        return ChannelStatus(
            "dispersion", BLIND, floor, "k_ex s^-1",
            f"a {reference_population:.0%} state needs k_ex >= {floor:.1e} s^-1 to be "
            f"fitted in {total_time_ns:.0f} ns, but cryptic opening runs near "
            f"{reference_k_ex:.0e} s^-1")
    if getattr(fit, "reliable", False):
        return ChannelStatus("dispersion", SEEN, floor, "k_ex s^-1")
    return ChannelStatus("dispersion", NO_SIGNAL, floor, "k_ex s^-1",
                         getattr(fit, "note", "") or "no reliable two-state fit")


def tail_status(x: np.ndarray, total_time_ns: float, *,
                reference_k_ex: float = REFERENCE_K_EX,
                reference_population: float = REFERENCE_POPULATION,
                min_clusters: int = 25) -> ChannelStatus:
    """Status of the extreme-value channel for one residue.

    A refused fit is always ``blind``, never ``no_signal``.  :func:`fit_tail`
    walks its threshold down to the median before giving up, so a refusal says
    the trajectory could not supply enough independent excursions -- a statement
    about the sampling, not about the site.

    A *successful* fit is not by itself evidence that anything opened: a
    generalised Pareto fits the tail of ordinary thermal noise perfectly well.
    Whether the return level is large is the score's business, not this
    function's.
    """
    floor = tail_floor(total_time_ns, reference_population, min_clusters)
    f = fit_tail(np.asarray(x, dtype=float).ravel())
    if f.ok:
        return ChannelStatus("tail", SEEN, floor, "k_ex s^-1")
    return ChannelStatus(
        "tail", BLIND, floor, "k_ex s^-1",
        f"{f.note or 'fit refused'}; the threshold was walked to the median and "
        f"still could not supply {min_clusters} independent excursions")


def breathing_status(contacts: np.ndarray, *,
                     reference_population: float = REFERENCE_POPULATION,
                     min_open_frames: int = 10) -> ChannelStatus:
    """Status of the protection-factor channel for one residue.

    The floor is a property of the trajectory's length alone, deliberately not
    of the residue's own autocorrelation time.  A residue that opens in long
    episodes has a long correlation time *because* of the opening, so charging
    it for that would make the observability test circular -- it would declare
    the strongest sites unobservable on the strength of their own signal.  That
    is exactly what an earlier version of this function did, and it sent the one
    genuinely observable cryptic site in the demo to last place.

    Being a time average with no fit to support, this channel needs occupied
    frames rather than completed round trips, which is why it keeps reach where
    the other two lose it.
    """
    x = np.asarray(contacts, dtype=float).ravel()
    floor = breathing_floor(x.size, 1.0, min_open_frames)
    if floor > reference_population:
        return ChannelStatus(
            "breathing", BLIND, floor, "population",
            f"{x.size} frames hold at most {floor:.3g} resolvable population, "
            f"above the {reference_population:.0%} reference")
    return ChannelStatus("breathing", SEEN, floor, "population")


# --------------------------------------------------------------------------
# Residue outcome and cohort accounting.
# --------------------------------------------------------------------------

def residue_outcome(channels: dict) -> str:
    """``SCORED`` if anything was seen, ``NEGATIVE`` if a sighted channel found
    nothing, ``ABSENT`` if every channel was blind."""
    if any(c.status == SEEN for c in channels.values()):
        return SCORED
    if any(c.status == NO_SIGNAL for c in channels.values()):
        return NEGATIVE
    return ABSENT


def assess(residue: int, *, fit=None, openness: Optional[np.ndarray] = None,
           contacts: Optional[np.ndarray] = None,
           total_time_ns: float = 100.0, **kw) -> ResidueStatus:
    """Assemble the per-channel verdicts for one residue."""
    channels: dict = {}
    if fit is not None:
        channels["dispersion"] = dispersion_status(fit, total_time_ns, **kw)
    if openness is not None:
        channels["tail"] = tail_status(openness, total_time_ns, **{
            k: v for k, v in kw.items() if k in
            ("reference_k_ex", "reference_population", "min_clusters")})
    if contacts is not None:
        channels["breathing"] = breathing_status(contacts, **{
            k: v for k, v in kw.items() if k in
            ("reference_population", "min_open_frames")})
    return ResidueStatus(residue, residue_outcome(channels), channels)


def summarise(statuses: Sequence[ResidueStatus]) -> dict:
    """Abstention accounting for a whole protein.

    A residue whose every channel was blind must be counted, not dropped.
    Reporting an enrichment over the observable subset alone states it
    conditional on observability, and observability here is correlated with the
    thing being ranked: the slower a site opens, the more likely it is both a
    genuine cryptic pocket and invisible.
    """
    n = len(statuses)
    counts = {SCORED: 0, NEGATIVE: 0, ABSENT: 0}
    blind: dict[str, int] = {}
    for s in statuses:
        counts[s.outcome] = counts.get(s.outcome, 0) + 1
        for name in s.blind_channels:
            blind[name] = blind.get(name, 0) + 1
    return {
        "n_residues": n,
        "counts": counts,
        "absent_share": (counts[ABSENT] / n) if n else 0.0,
        "blind_by_channel": blind,
    }
