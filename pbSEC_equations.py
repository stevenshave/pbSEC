from mpmath import mpf, sqrt, power, mp

mp.dps = 100


def one_to_one_binding(p: float, l: float, kdpl: float):
    """Return protein-ligand complex amount as a float

    Args:
        p (float): Protein concentration
        l (float): Ligand concentration
        kdpl (float): Protein-ligand KD

    Returns:
        float: Complex concentration
    """
    p = mpf(p)
    l = mpf(l)
    kdpl = mpf(kdpl)
    return float((p + kdpl + l - sqrt(-4 * p * l + power(p + kdpl + l, 2))) / 2.0)


def one_to_one_binding_mpf(p: mpf, l: mpf, kdpl: mpf):
    """Return protein-ligand complex amount as an mpf

    Args:
        p (mpf): Protein concentration
        l (mpf): Ligand concentration
        kdpl (mpf): Protein-ligand KD

    Returns:
        mpf: Complex concentration
    """
    return (p + kdpl + l - sqrt(-4 * p * l + power(p + kdpl + l, 2))) / 2.0


def pbSEC(p: float, l: float, kdpl: float, recovery_efficiency: float = 1.0):
    """Calculate the amount protein-ligand complex recovered in a pbSEC run

    Args:
        p (float): Protein concentration
        l (float): Ligand concentration
        kdpl (float): Protein-ligand dissociation constant (KD)
        recovery_efficiency (float, optional): Fraction protein recovered on each iteration. Defaults to 1.0.
    Returns:
        float: amount complex

    """
    p, l, kdpl = mpf(p), mpf(l), mpf(kdpl)
    pl_complex = one_to_one_binding_mpf(p, l, kdpl) * recovery_efficiency
    return float(pl_complex)


def pbSEC_simulate_n_rounds(
    p: float, l: float, kdpl: float, n: int, recovery_efficiency: float = 1.0
):
    """Calculate the amount protein-ligand complex recovered in iterative pbSEC runs and return as a list of floats

    Args:
        p (float): Protein concentration
        l (float): Ligand concentration
        kdpl (float): Protein-ligand dissociation constant (KD)
        n (int): Number of pbSEC rounds
        recovery_efficiency (float, optional): Fraction protein recovered on each iteration. Defaults to 1.0.
    Returns:
        list(float): List of complex amounts at each iteration

    """
    p, l, kdpl = mpf(p), mpf(l), mpf(kdpl)
    round_count = 0
    complex_concentrations = []
    for _ in range(n):
        l = one_to_one_binding_mpf(p, l, kdpl) * recovery_efficiency
        p = p * recovery_efficiency
        complex_concentrations.append(float(l))
    return complex_concentrations


def pbSEC_iterate_until_undetectable(
    p: float,
    l: float,
    kdpl: float,
    l_detection_limit: float,
    recovery_efficiency: float = 1.0,
):
    """Calculate amount returned in interative pbSEC rounds until a certain detection limit is reached.

    Args:
        p (float): Protein concentration.
        l (float): Ligand concentration.
        kdpl (float): Protein-ligand dissociation constant (KD).
        l_detection_limit (float): Ligand detection limit.
        recovery_efficiency (float, optional): Fraction protein recovered on each iteration. Defaults to 1.0.
    Returns:
        list(float): List of complex amounts at each iteration until undetectable.

    """
    p, l, kdpl = mpf(p), mpf(l), mpf(kdpl)
    round_count = 0
    complex_concentrations = []
    while l >= l_detection_limit:
        l = one_to_one_binding_mpf(p, l, kdpl) * recovery_efficiency
        p = p * recovery_efficiency
        complex_concentrations.append(float(l))
    return complex_concentrations[:-1]
