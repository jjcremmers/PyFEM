# pyfem/api.py

'''
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Union

# If you prefer, validate with pydantic (optional)
# from pydantic import BaseModel

Props = Dict[str, Any]  # or a TypedDict / Pydantic model
Output = Dict[str, Any] # or a Results dataclass

def _load_props(props: Union[str, Path, Props]) -> Props:
    """Accepts a dict, JSON, or YAML file path."""
    if isinstance(props, (str, Path)):
        p = Path(props)
        if p.suffix.lower() in {".yml", ".yaml"}:
            import yaml
            return yaml.safe_load(p.read_text())
        elif p.suffix.lower() == ".json":
            import json
            return json.loads(p.read_text())
        else:
            raise ValueError(f"Unsupported props file: {p}")
    elif isinstance(props, dict):
        return props
    else:
        raise TypeError("props must be a dict or a path to a YAML/JSON file")
'''

def run( a ):
    """Run a PyFEM analysis from a dict or YAML/JSON file path.

    Parameters
    ----------
    props : dict | str | Path
        Analysis settings (or a path to YAML/JSON).

    Returns
    -------
    dict
        A results dictionary (or switch to a dataclass).
    """

    print("A = ",a)
    
    return 2*a

