from .synthetic_data import (generate_synthetic_data,
                             generate_background,
                             generate_signals,
                            )


def test_generate_shape():
    data = generate_synthetic_data(shape=(10, 10))
    assert data.shape == (10, 10)


def test_background_level():
    background_level = 100
    background_sigma = 0.01
    data = generate_background(
        shape=(10, 10),
        background_level=background_level,
        background_sigma=background_sigma,
    )

    expected_min = 95
    expected_max = 105
    assert expected_min < data[0][0] < expected_max


def test_generate_deltas():
    positions = [(4, 4)]
    data = generate_signals(positions=positions, shape=(10, 10))
    assert data[positions[0]] > 0.
    assert data[positions[0][1], positions[0][0] - 1] == 0.
    assert data[0, 0] == 0.

def test_generate_multiple_deltas():
    positions = [(4, 4), (8, 8), (5, 5)]
    data = generate_signals(positions=positions, shape=(10, 10))
    assert data[0, 0] == 0.
    for position in positions:
        assert data[position] > 0.
        assert data[position[1], position[0] - 1] == 0.
