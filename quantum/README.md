## Installation and Usage with `uv`

[`uv`](https://github.com/astral-sh/uv) is a fast Python package and environment manager.

### Requirements

* `uv` installed (choose one):

  ```bash
  # Recommended installer
  curl -LsSf https://astral.sh/uv/install.sh | sh
  ```
  Or install via **pip**:

  ```bash
  pip install uv
  ```

### Setup

```bash
    cd quantum-protein-folding/src/quantum
    uv sync
```

This creates a virtual environment and installs dependencies from
`pyproject.toml`

### Running

```
uv run main.py
```

### Adding Dependencies

```
uv add <package-name>
```


### Running tests

```
cd quantum
```

```
uv run pytest -v -s
```