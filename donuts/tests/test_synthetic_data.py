from .synthetic_data import generate_synthetic_data


def test_generate_shape():
    data = generate_synthetic_data(shape=(10, 10))
    assert data.shape == (10, 10)
