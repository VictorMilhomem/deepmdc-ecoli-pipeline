import os
import json
import pandas as pd
from functools import partial
from os import path, listdir, mkdir
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, Optional

from src.utils.handle_files import write_table, read_tsv
from pipeline_regulatory.processing.handle_regulatory_processing import (
    process_json_regulatory,
    process_tsv_regulatory,
    read_json_any
)


def create_regulatory_elements_table(
    metadata_file: str = "database/metadata.tsv",
    complete_table_dir: str = "database/complete_tables_regulatory",
    jsons_dir: str = "jsons",
    regulatory_tables_dir: str = "database/regulatory_elements",
    catalog_json_output: str = "database/all_elements.json",
    sequence_mode: str = "nt",
    include_srnas: bool = True,
    filter_orfs: bool = False
):
    """
    Gera tabelas de elementos genômicos regulatórios a partir dos JSONs do Bakta.

    Diferenças em relação ao 'create_orf_table' original:
    - Processa ORFs e sRNAs (não descarta sequências não-codificantes).
    - Gera tabelas prontas para o DeepRC em 'regulatory_tables_dir'.
    - Cria o catálogo global 'all_elements.json' consolidado com informações
      biológicas de cada sequência para a análise de pesos de atenção posterior.

    Args:
        metadata_file (str): TSV com os IDs das amostras a processar.
        complete_table_dir (str): Pasta onde salvar tabelas intermediárias detalhadas.
        jsons_dir (str): Pasta contendo os JSONs (.json ou .json.xz) do Bakta.
        regulatory_tables_dir (str): Pasta de saída final dos TSVs (orf / templates).
        catalog_json_output (str): Caminho para salvar o catálogo global desduplicado.
        sequence_mode (str): 'nt' (DNA, recomendado para incluir sRNAs/promotores) ou 'aa'.
        include_srnas (bool): Se True, inclui pequenos RNAs regulatórios anotados.
        filter_orfs (bool): Se True, descarta proteínas hipotéticas.
    """
    metadata_file = path.abspath(metadata_file)
    metadata = pd.read_csv(metadata_file, sep='\t')
    ids_to_use = {sample_id.split(".")[0] for sample_id in metadata['ID']}

    complete_table_dir = path.abspath(complete_table_dir)
    regulatory_tables_dir = path.abspath(regulatory_tables_dir)
    if not path.exists(complete_table_dir):
        os.makedirs(complete_table_dir, exist_ok=True)
    if not path.exists(regulatory_tables_dir):
        os.makedirs(regulatory_tables_dir, exist_ok=True)

    base_dir = path.abspath(jsons_dir)
    if not path.exists(base_dir):
        print(f"[Aviso] Pasta {base_dir} não encontrada. Certifique-se de disponibilizar os JSONs.")
        return

    # Localizar arquivos compatíveis (.json ou .json.xz)
    json_files = []
    for f in listdir(base_dir):
        clean_id = f.replace(".json.xz", "").replace(".json", "")
        if clean_id in ids_to_use:
            json_files.append(path.join(base_dir, f))

    if len(json_files) == 0:
        print(f"[Aviso] Nenhum arquivo JSON correspondente aos IDs do metadata foi encontrado em {base_dir}.")
        return

    print(f"[Regulatory] Processando {len(json_files)} arquivos genômicos (Modo: {sequence_mode}, sRNAs: {include_srnas})...")

    worker = partial(
        json_to_regulatory_table,
        complete_table_dir=complete_table_dir,
        regulatory_tables_dir=regulatory_tables_dir,
        sequence_mode=sequence_mode,
        include_srnas=include_srnas,
        filter_orfs=filter_orfs
    )

    with ProcessPoolExecutor() as executor:
        executor.map(worker, json_files)

    print(f"[Regulatory] Tabelas de elementos geradas em: {regulatory_tables_dir}")

    # Gerar catálogo biológico consolidado (all_elements.json)
    if catalog_json_output:
        print(f"[Regulatory] Consolidando catálogo biológico em: {catalog_json_output}...")
        build_catalog_json(
            complete_table_dir=complete_table_dir,
            output_catalog_path=catalog_json_output
        )


def json_to_regulatory_table(
    json_file: str,
    complete_table_dir: str,
    regulatory_tables_dir: str,
    sequence_mode: str = "nt",
    include_srnas: bool = True,
    filter_orfs: bool = False
):
    """
    Processa um único JSON do Bakta, gerando a tabela completa e o TSV final.
    """
    rows = process_json_regulatory(
        json_file=json_file,
        sequence_mode=sequence_mode,
        include_srnas=include_srnas,
        filter_orfs=filter_orfs
    )

    base_name = path.basename(json_file).replace(".json.xz", "").replace(".json", "")
    table_name = f"{base_name}.tsv"
    complete_table = path.join(complete_table_dir, table_name)

    write_table(complete_table, rows)

    reg_table = path.abspath(path.join(regulatory_tables_dir, table_name))
    df = read_tsv(complete_table)

    processed_df = process_tsv_regulatory(df)
    write_table(reg_table, processed_df)


def build_catalog_json(complete_table_dir: str, output_catalog_path: str):
    """
    [NOVO] Cria o catálogo biológico consolidado (all_elements.json) a partir
    das tabelas completas geradas. Armazena:
    {
        "<sequencia>": {
            "type": "cds" | "srna" | "promoter",
            "product": "...",
            "locus": "..."
        }
    }
    """
    catalog: Dict[str, dict] = {}
    tsv_files = [path.join(complete_table_dir, f) for f in listdir(complete_table_dir) if f.endswith(".tsv")]

    for tsv_f in tsv_files:
        df = pd.read_csv(tsv_f, sep="\t")
        for _, row in df.iterrows():
            seq = row["sequence"]
            if seq not in catalog:
                catalog[seq] = {
                    "type": row.get("type", "unknown"),
                    "product": row.get("product", ""),
                    "locus": row.get("locus", "")
                }

    os.makedirs(path.dirname(output_catalog_path), exist_ok=True)
    with open(output_catalog_path, "w", encoding="utf-8") as out:
        json.dump(catalog, out, indent=2)

    print(f"[Regulatory] Catálogo biológico salvo com {len(catalog)} sequências únicas.")


if __name__ == "__main__":
    create_regulatory_elements_table()

