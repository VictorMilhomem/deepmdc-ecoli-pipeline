import os
import lzma
import json
import pandas as pd
from datetime import datetime
from collections import Counter
from typing import List, Union, Callable, Tuple
import src.utils.handle_files as hf
from torch.utils.data import DataLoader
from deeprc.task_definitions import TaskDefinition
import h5py
from deeprc.dataset_converters import DatasetToHDF5
from deeprc.dataset_readers import RepertoireDataset, no_stack_collate_fn, \
    no_sequence_count_scaling

# ==============================================================================
# PIPELINE REGULATÓRIO - PROCESSAMENTO DE DADOS
# Adaptado de: src/utils/handle_processing.py
# Principais mudanças:
# 1. Suporte à extração conjunta de ORFs + regiões regulatórias (sRNAs como MicF,
#    promotores e regiões upstream).
# 2. Suporte à representação em nucleotídeos (DNA: A, C, G, T) ou aminoácidos (AA).
# 3. Leitura direta de JSONs (.json ou .json.xz do Bakta).
# 4. keep_dataset_in_ram padronizado como False (corrige erro de 33.5 GB de RAM).
# 5. Parametrização do alfabeto de caracteres (sequence_characters) no DatasetToHDF5.
# ==============================================================================

DNA_NUCLEOTIDES = ('A', 'C', 'G', 'T')
AMINO_ACIDS = (
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
)


def read_json_any(json_path: str) -> dict:
    """
    Lê arquivos JSON tanto em formato simples (.json) quanto
    comprimido via LZMA (.json.xz).
    """
    if not os.path.exists(json_path):
        return {}

    if json_path.endswith(".xz"):
        with lzma.open(json_path, "rt", encoding="utf-8") as inp:
            return json.load(inp)
    else:
        with open(json_path, "r", encoding="utf-8") as inp:
            return json.load(inp)


def process_json_regulatory(
    json_file: str,
    sequence_mode: str = "nt",
    include_srnas: bool = True,
    filter_orfs: bool = False,
    min_len: int = 15,
    max_len: int = 5000
) -> List[List[Union[str, int]]]:
    """
    Processa um JSON do Bakta extraindo tanto ORFs codificadoras quanto
    regiões regulatórias (sRNAs, ncRNAs) e sequências não codificantes.

    Diferenças em relação ao original:
    - No original: extraía apenas feature.get("aa", "") de ORFs proteicas.
    - No regulatório:
        * Se sequence_mode == "nt": extrai nucleotídeos (DNA) tanto de CDS
          quanto de sRNAs/ncRNAs, unificando o alfabeto em 4 letras (A, C, G, T).
        * Se sequence_mode == "aa": extrai aminoácidos para CDS e nucleotídeos
          para sRNAs (modo híbrido).
        * Identifica sRNAs e ncRNAs (ex: MicF, CsrB, RyhB, etc.) presentes
          na chave 'features' do Bakta.

    Args:
        json_file (str): Caminho para o JSON (.json ou .json.xz) gerado pelo Bakta.
        sequence_mode (str): 'nt' para nucleotídeos (recomendado para integrar
                             promotores e sRNAs) ou 'aa' para aminoácidos.
        include_srnas (bool): Se True, inclui pequenos RNAs regulatórios anotados.
        filter_orfs (bool): Se True, filtra ORFs com anotação hipotética.
        min_len (int): Tamanho mínimo da sequência a ser incluída.
        max_len (int): Tamanho máximo da sequência a ser incluída.

    Returns:
        List[List[Union[str, int]]]: Linhas contendo [locus, type, product, sequence, count].
    """
    data = read_json_any(json_file)
    features = data.get("features", [])

    sequences = []
    rows = []

    SRNA_TYPES = {"srna", "ncrna", "rrna", "trna", "crispr"}

    for feature in features:
        locus = feature.get("locus", feature.get("id", ""))
        feat_type = str(feature.get("type", "")).lower()
        product = feature.get("product", "")
        gene = feature.get("gene", "")

        # 1. Processar CDS (ORFs codificadoras de proteínas)
        if feat_type == "cds":
            if sequence_mode == "nt":
                seq = feature.get("nt", "").upper()
            else:
                seq = feature.get("aa", "").upper()

            if not seq:
                continue

            if filter_orfs and product == "hypothetical protein":
                continue

            if min_len <= len(seq) <= max_len:
                sequences.append(seq)
                rows.append([locus, "cds", product if product else gene, seq])

        # 2. Processar sRNAs e outros RNAs regulatórios não-codificantes
        elif include_srnas and feat_type in SRNA_TYPES:
            seq = feature.get("nt", "").upper().replace("U", "T")

            if not seq:
                continue

            if min_len <= len(seq) <= max_len:
                label_product = f"sRNA:{gene}" if gene else f"sRNA:{product}"
                sequences.append(seq)
                rows.append([locus, "srna", label_product, seq])

    # Contagem de ocorrências de cada sequência única no genoma (templates)
    seq_counter = Counter(sequences)

    for row in rows:
        row.append(seq_counter.get(row[3], 1))

    return rows


def process_tsv_regulatory(df: pd.DataFrame) -> pd.DataFrame:
    """
    Processa o DataFrame intermediário mantendo as colunas 'orf' e 'templates'
    necessárias para a entrada do DeepRC.
    """
    columns = {"sequence": "orf", "count": "templates"}
    selected_df = df[["sequence", "count"]].copy()
    selected_df.rename(columns=columns, inplace=True)
    return selected_df


def make_dataloader_regulatory(
    task_definition: TaskDefinition,
    metadata_file: str,
    genomesdata_path: str,
    n_worker_processes: int = 4,
    batch_size: int = 1,
    inputformat: str = 'NCL',
    keep_dataset_in_ram: bool = False,  
    ask_for_input: bool = False,
    metadata_file_id_column: str = 'ID',
    metadata_file_column_sep: str = '\t',
    sequence_column: str = "orf",
    sequence_counts_column: str = "templates",
    genome_files_column_sep: str = '\t',
    filename_extension: str = '.tsv',
    h5py_dict: dict = None,
    sequence_characters: Tuple[str, ...] = DNA_NUCLEOTIDES,  
    sequence_counts_scaling_fn: Callable = no_sequence_count_scaling,
    verbose: bool = True
) -> DataLoader:
    """
    Cria DataLoader para avaliação/inferência de elementos genômicos regulatórios.

    Diferenças principais:
    - keep_dataset_in_ram padrão é FALSE (leitura sob demanda em disco / lazy loading).
    - sequence_characters configurável (permite DNA com 'A','C','G','T' ou 20 AAs).
    """
    try:
        with h5py.File(genomesdata_path, 'r') as hf:
            n_genomes = hf['metadata']['n_samples'][()]
        hdf5_file = genomesdata_path
    except Exception:
        hdf5_file = genomesdata_path + ".hdf5"
        if verbose:
            print(f"[Regulatory] Convertendo: {genomesdata_path} -> {hdf5_file}")
            print(f"[Regulatory] Alfabeto utilizado: {sequence_characters}")

        converter = DatasetToHDF5(
            repertoiresdata_directory=genomesdata_path,
            sequence_column=sequence_column,
            sequence_counts_column=sequence_counts_column,
            column_sep=genome_files_column_sep,
            filename_extension=filename_extension,
            sequence_characters=sequence_characters,
            h5py_dict=h5py_dict,
            verbose=verbose
        )
        converter.save_data_to_file(
            output_file=hdf5_file, n_workers=n_worker_processes
        )
        with h5py.File(hdf5_file, 'r') as hf:
            n_genomes = hf['metadata']['n_samples'][()]

    dataset = RepertoireDataset(
        metadata_filepath=metadata_file,
        hdf5_filepath=hdf5_file,
        sample_id_column=metadata_file_id_column,
        metadata_file_column_sep=metadata_file_column_sep,
        task_definition=task_definition,
        keep_in_ram=keep_dataset_in_ram,
        inputformat=inputformat,
        sequence_counts_scaling_fn=sequence_counts_scaling_fn
    )

    dataloader = DataLoader(
        dataset=dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=1,
        collate_fn=no_stack_collate_fn
    )
    return dataloader

