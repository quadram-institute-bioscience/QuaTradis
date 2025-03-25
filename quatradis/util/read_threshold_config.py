import yaml

def load_thresholds(config_path):
    """Load YAML thresholds dynamically from a given path."""
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)
    return config