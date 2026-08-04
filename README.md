## DeepMDC-RawData

This repository contains the code to download the raw data used in the [DeepMDC](https://github.com/CABGenOrg/deepmdc) project.


### Prerequisites

- Conda or Miniconda/Anaconda installed
- Python available through Conda
- Poetry installed

### Using Conda

Create and activate a new environment:

```bash
conda create -n deepmdc-rawdata python=3.11
conda activate deepmdc-rawdata
```

Install the Dependencies listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

### Using Poetry

Install the project dependencies:

```bash
poetry install
```


### Run the project

Run the application through **Poetry**:

```bash
poetry run python src/deepmdc_rawdata/get_data.py
```

Run the application through **Conda**:

```bash
python src/deepmdc_rawdata/get_data.py
```

