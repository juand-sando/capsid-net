import json
from pathlib import Path


def write_preprocess_config(path, payload):
    output_path = Path(path)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def load_preprocess_config(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
