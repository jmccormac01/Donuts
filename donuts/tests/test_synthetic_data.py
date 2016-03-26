from .synthetic_data import generate_synthetic_data


def test_generate_shape():
    data = generate_synthetic_data(shape=(10, 10))
    assert data.shape == (10, 10)


def test_background_level():
    background_level = 100
    background_sigma = 0.01
    data = generate_synthetic_data(
        background_level=background_level,
        background_sigma=background_sigma,
    )

    expected_min = 95
    expected_max = 105
    assert expected_min < data[0][0] < expected_max

