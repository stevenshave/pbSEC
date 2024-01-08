"""
Functions to simulate pbSEC rounds

"""

__all__ = [
    'one_to_one_binding',
    'pbSEC',
    'pbSEC_simulate_n_rounds',
    'pbSEC_iterate_until_undetectable',
]

from mpmath import mpf, sqrt, power, mp

mp.dps = 100


def one_to_one_binding(p0: float, l0: float, kdpl: float) -> float:
    """Get protein-ligand complex concentration as a float

    Parameters
    ----------
    p0 : float
        Initial free protein concentration
    l0 : float
        Initial free ligand concentration
    kdpl : float
        Protein ligand affinity expressed as a dissociation constant (KD)

    Returns
    -------
    float
        Protein-ligand complex concentration
    """
    p0 = mpf(p0)
    l0 = mpf(l0)
    kdpl = mpf(kdpl)
    return float((p0 + kdpl + l0 - sqrt(-4 * p0 * l0 + power(p0 + kdpl + l0, 2))) / 2.0)


def _one_to_one_binding_mpf(p0: mpf, l0: mpf, kdpl: mpf) -> mpf:
    """Get protein-ligand complex concentration as a mpmath high precision float

    Parameters
    ----------
    p0 : float
        Initial free protein concentration
    l0 : float
        Initial free ligand concentration
    kdpl : float
        Protein ligand affinity expressed as a dissociation constant (KD)

    Returns
    -------
    mpmath.mpf
        Protein-ligand complex concentration
    """
    return (p0 + kdpl + l0 - sqrt(-4 * p0 * l0 + power(p0 + kdpl + l0, 2))) / 2.0


def pbSEC(p0: float, l0: float, kdpl: float, recovery_efficiency: float = 1.0) -> float:
    """Calculate complex concentration recovered recovered from a pbSEC run

    Parameters
    ----------
    p0 : float
        Initial free protein concentration
    l0 : float
        Initial free ligand concentration
    kdpl : float
        Protein ligand affinity expressed as a dissociation constant (KD)
    recovery_efficiency : float, optional
        Fraction of protein complex recovered after each pbSEC iteration. 100 %
        recovery should be expressed as 1.0, similarly 80 % recovery as 0.8. By
        default 1.0

    Returns
    -------
    float
        Concentration of recovered protein-ligand complex
    """
    p0, l0, kdpl = mpf(p0), mpf(l0), mpf(kdpl)
    pl_complex = _one_to_one_binding_mpf(p0, l0, kdpl) * recovery_efficiency
    return float(pl_complex)


def pbSEC_simulate_n_rounds(
    p0: float,
    l0: float,
    kdpl: float,
    num_iterations: int,
    recovery_efficiency: float = 1.0,
) -> float:
    """Calc complex concentration recovered after iterative pbSEC runs

    Parameters
    ----------
    p0 : float
        Initial free protein concentration
    l0 : float
        Initial free ligand concentration
    kdpl : float
        Protein ligand affinity expressed as a dissociation constant (KD)
    num_iterations : int
        _description_
    recovery_efficiency : float, optional
        Fraction of protein complex recovered after each pbSEC iteration. 100 %
        recovery should be expressed as 1.0, similarly 80 % recovery as 0.8. By
        default 1.0

    Returns
    -------
    float
        Concentration of recovered protein-ligand complex after iterative rounds
    """

    p0, l0, kdpl = mpf(p0), mpf(l0), mpf(kdpl)
    round_count = 0
    complex_concentrations = []
    for _ in range(num_iterations):
        l0 = _one_to_one_binding_mpf(p0, l0, kdpl) * recovery_efficiency
        p0 = p0 * recovery_efficiency
        complex_concentrations.append(float(l0))
    return complex_concentrations


def pbSEC_iterate_until_undetectable(
    p0: float,
    l0: float,
    kdpl: float,
    l_detection_limit: float,
    recovery_efficiency: float = 1.0,
):
    """Iterate pbSEC rounds until ligand is undetectable

    Parameters
    ----------
    p0 : float
        Initial free protein concentration
    l0 : float
        Initial free ligand concentration
    kdpl : float
        Protein ligand affinity expressed as a dissociation constant (KD)
    l_detection_limit : float
        ligand detection limit
    recovery_efficiency : float, optional
        Fraction of protein complex recovered after each pbSEC iteration. 100 %
        recovery should be expressed as 1.0, similarly 80 % recovery as 0.8. By
        default 1.0

    Returns
    -------
    List[float, ...]
        List containing the complex concentration recovered at each pbSEC
        iteration until undetectable
    """
    p0, l0, kdpl = mpf(p0), mpf(l0), mpf(kdpl)
    round_count = 0
    complex_concentrations = []
    while l0 >= l_detection_limit:
        l0 = _one_to_one_binding_mpf(p0, l0, kdpl) * recovery_efficiency
        p0 = p0 * recovery_efficiency
        complex_concentrations.append(float(l0))
    return complex_concentrations[:-1]
