import numpy as np
from .helpers import (
    Bunch,
    fort_write,
    find_indices_with_tol,
)
from .constants import MAX_NUM_POINTS
from .datamodels.models import Dataset


def propagate_prior_to_dataset(dataset: Dataset, reaction_prior):
    dataset = Bunch(dataset.dict(use_arrays=True))
    reaction_prior = reaction_prior.dict()
    prior_energy_mesh = [p['energies'] for p in reaction_prior.values()]
    prior_cross_section = [p['cross_sections'] for p in reaction_prior.values()]

    if dataset.quantity_type in (1, 2):
        E11, E22, EQ, Q, mxm1 = propagate_prior_to_CS_and_CS_SHAPE_dataset(
            dataset, prior_energy_mesh, prior_cross_section
        )
    elif dataset.quantity_type in (3, 4):
        E11, E22, EQ, Q, mxm1 = propagate_prior_to_RATIO_and_RATIO_SHAPE_dataset(
            dataset,  prior_energy_mesh, prior_cross_section
        )
    elif dataset.quantity_type in (5, 8):
        E11, E22, EQ, Q, mxm1 = propagate_prior_to_SUM_and_SHAPE_OF_SUM_dataset(
            dataset, prior_energy_mesh, prior_cross_section
        )
    elif dataset.quantity_type in (7, 9):
        E11, E22, EQ, Q, mxm1 = propagate_prior_to_CS_VS_SUM_PLUS_SHAPE_dataset(
            dataset, prior_energy_mesh, prior_cross_section
        )
    elif dataset.quantity_type in (6,):
        E11, E22, EQ, Q, mxm1 = (None,) * 5

    else:
        raise ValueError(f'Invalid `quantity_type={dataset.quantity_type}`')
    assert dataset.quantity_type >= 1 and dataset.quantity_type <= 9

    return E11, E22, EQ, Q, mxm1


def propagate_prior_to_CS_and_CS_SHAPE_dataset(
    dataset, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((MAX_NUM_POINTS,), dtype=float)
    # CS + CS SHAPE
    M1 = dataset.reaction_ids[0] - 1
    NON = len(prior_energy_mesh[M1])
    NON1 = NON-1
    E11 = (prior_energy_mesh[M1][0] + prior_energy_mesh[M1][1]) / 2.
    E22 = (prior_energy_mesh[M1][NON-1] + prior_energy_mesh[M1][NON1-1]) / 2.
    for K in range(NON):
        EQ[K] = prior_energy_mesh[M1][K]
        Q[K] = prior_cross_section[M1][K]
    format4611 = "(i6,'  cr. sec. apriori for interp. ')"
    fort_write(None, format4611, [NON-1])
    return E11, E22, EQ, Q, NON-1


def propagate_prior_to_RATIO_and_RATIO_SHAPE_dataset(
    dataset, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((MAX_NUM_POINTS,), dtype=float)
    # RATIO + RATIO SHAPE
    M1 = dataset.reaction_ids[0] - 1
    M2 = dataset.reaction_ids[1] - 1
    NO1 = len(prior_energy_mesh[M1])
    NON = len(prior_energy_mesh[M2])
    mxm = 0

    for K in range(NO1):
        # find matching energies
        L = find_indices_with_tol(prior_energy_mesh[M2][:NON], prior_energy_mesh[M1][K:K+1], atol=1e-8, rtol=1e-4).item()
        if L == -1:  # not found
            continue
        EQ[mxm] = prior_energy_mesh[M1][K]
        if prior_cross_section[M2][L] == 0.:
            format4618 = ":(' apriori constr. ',4i5,f10.7,'*****************')"
            fort_write(None, format4618, [M1-1, K+1, M2-1, L+1, EQ[mxm]])
            exit()

        Q[mxm] = prior_cross_section[M1][K] / prior_cross_section[M2][L]
        mxm += 1

    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[mxm-1] + EQ[mxm-2]) / 2.

    format4612 = "(i6,'  ratio apriori for interp. ')"
    fort_write(None, format4612, [mxm])
    return E11, E22, EQ, Q, mxm-1


def propagate_prior_to_SUM_and_SHAPE_OF_SUM_dataset(
    dataset, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((MAX_NUM_POINTS,), dtype=float)
    # SUM AND SHAPE OF SUM
    M1 = dataset.reaction_ids[0] - 1
    M2 = dataset.reaction_ids[1] - 1
    M3 = dataset.reaction_ids[2] - 1
    NO1 = len(prior_energy_mesh[M1])
    NON = len(prior_energy_mesh[M2])
    if M3 != -1:
        NO3 = len(prior_energy_mesh[M3])

    mxm = 0
    for K in range(NO1):  # 23

        L = find_indices_with_tol(
            prior_energy_mesh[M2][:NON], prior_energy_mesh[M1][K:K+1], atol=1e-8, rtol=1e-4
        ).item()
        if L == -1:  # not found
            continue

        if M3 != -1:
            J = find_indices_with_tol(
                prior_energy_mesh[M3][:NO3], prior_energy_mesh[M1][K:K+1], atol=1e-8, rtol=1e-4
            ).item()
            if J == -1:  # not found
                continue

        EQ[mxm] = prior_energy_mesh[M1][K]
        Q[mxm] = prior_cross_section[M1][K] + prior_cross_section[M2][L]
        if M3 != -1:
            Q[mxm] = Q[mxm] + prior_cross_section[M3][J]

        mxm += 1

    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[mxm-1] + EQ[mxm-2]) / 2.

    format4613 = "(i6,'  sum apriori for interp. ')"
    fort_write(None, format4613, [mxm])
    return E11, E22, EQ, Q, mxm-1


def propagate_prior_to_CS_VS_SUM_PLUS_SHAPE_dataset(
    dataset, prior_energy_mesh, prior_cross_section
):
    EQ = np.empty((200,), dtype=float)
    Q = np.zeros((MAX_NUM_POINTS,), dtype=float)
    #  RATIO OF CS VS. SUM + SHAPE
    M1 = dataset.reaction_ids[0] - 1
    M2 = dataset.reaction_ids[1] - 1
    M3 = dataset.reaction_ids[2] - 1
    NO1 = len(prior_energy_mesh[M1])
    NON = len(prior_energy_mesh[M2])
    NO3 = len(prior_energy_mesh[M3])
    mxm = 0
    pymxm = mxm - 1

    for K in range(NO1):

        L = find_indices_with_tol(
            prior_energy_mesh[M2][:NON], prior_energy_mesh[M1][K:K+1], atol=1e-8, rtol=1e-4
        ).item()
        if L == -1:  # not found
            continue

        if M3 != -1:
            J = find_indices_with_tol(
                prior_energy_mesh[M3][:NO3], prior_energy_mesh[M1][K:K+1], atol=1e-8, rtol=1e-4
            ).item()
            if J == -1:  # not found
                continue

        EQ[mxm] = prior_energy_mesh[M1][K]
        Q[mxm] = (
            prior_cross_section[M1][K] /
            (prior_cross_section[M2][L] + prior_cross_section[M3][J])
        )
        mxm += 1

    mxm1 = mxm - 1
    E11 = (EQ[0] + EQ[1]) / 2.
    E22 = (EQ[pymxm] + EQ[mxm1-1]) / 2.

    format4614 = "(i6,'  cr. sec. vs sum apriori for interp. ')"
    fort_write(None, format4614, [mxm])
    return E11, E22, EQ, Q, mxm1
